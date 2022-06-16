%% Testing 
[stackorg, pathname, index] = uigetfile('.tif', 'Open CELL Image');
stacktiff = tiffreadVolume(stackorg);

stacktiff_thrcut = stacktiff(:,20:400,20:100, :);

parameters1 = [50,500,0.2,0.5]; %volumethreshmin1,volumethreshmax1, extentthreshmin1,extentthreshmax1
parameters2 = [50]; %volumethreshmin2

threshstack(:,:,:,1) = imbinarize(stacktiff_thrcut(:,:,:,1)); %nuclei
threshstack(:,:,:,2) = imbinarize(stacktiff_thrcut(:,:,:,2)); %actin


%[nuccount, actincount, S1totalsort, labeledstack] = sprout_analysis(stacktiff_thrcut, threshstack,parameters1, parameters2);

%%

%function [nuccount, actincount, S1totalsort, labeledstack] = nuclei_analysis(stacktiff_thrcut, threshstack,parameters1, parameters2, actin_mean)
    
%parameters for threshholding
    volumethreshmin1 = parameters1(1); 
    volumethreshmax1 = parameters1(2);
    extentthreshmin1 = parameters1(3);
    extentthreshmax1 = parameters1(4);
    
    % Hole Filling

    threshstack(:,:,:,1) = imfill(threshstack(:,:,:,1), 'holes');
    
    %Pre-analysis

    CC1 = bwconncomp(threshstack(:,:,:,1), 26);



    S1idlist = regionprops3(CC1,'VoxelList');
    S1vol = regionprops3(CC1, 'Volume', 'Extent'); %Generating table based on parameters
    tableS1 = table2array(S1vol);


    filtCC1 = false(size(threshstack(:,:,:,1)));
    nuccount = 0;
    
    %% Nuclei loop
    for i=1:(size(S1vol,1))
        if tableS1(i,1) >= volumethreshmin1 && tableS1(i,1) <= volumethreshmax1 && tableS1(i,2) >= extentthreshmin1 && tableS1(i,2) <= extentthreshmax1
            specind = cell2mat(S1idlist{i,1});
            y = specind(:,1);
            x = specind(:,2);
            z = specind(:,3);
            for j = 1:(size(specind,1))
              filtCC1(x(j),y(j),z(j)) = true;
            end            
            nuccount = nuccount + 1;
            clear specind;
        end
    end
   
 
    %% New tables and initialization of nucleiRGBstack
    S1total = regionprops3(filtCC1,'Volume','Centroid','Orientation','PrincipalAxisLength','SurfaceArea','Extent','Orientation', 'EigenValues', 'EigenVectors',"EquivDiameter");
    S2total = regionprops3(filtCC2,'Volume','Centroid','Orientation','PrincipalAxisLength','SurfaceArea','Extent','Orientation', 'EigenValues', 'EigenVectors',"EquivDiameter");
    
    S1total = splitvars(S1total,2); %sorting rows by z value
    S1sorted = sortrows(S1total,4);
    
    namelist1 = linspace(1,size(S1total,1), size(S1total,1));
    namelist1 = permute(namelist1, [2,1]);
    tablelist1 = table(namelist1);
    tablelist1.Properties.VariableNames(1) = "Number";
    S1totalsort = [tablelist1 S1sorted];

    nucleiRGBstack(:,:,:,1) = filtCC1 * 255; %binarized layer
    
    nucleiRGBstack(:,:,:,2) = im2double(stacktiff_thrcut(:,:,:,2));
    
    nucleiblue = zeros(size(filtCC1));
    
    %% Identifying nuclei, setting them to 255 in blue stack
    nucloctable = [];
    for i = 1:nuccount
        nucnum = S1totalsort{i,1};
        x = round(S1totalsort{i,3});
        y = round(S1totalsort{i,4});
        z = round(S1totalsort{i,5});


        if x < 1
            x = 1;
        end
        if y < 1
            y = 1;
        end
        if z < 1
            z = 1;
        end
        if x > size(filtCC1,2)
            x = size(filtCC1,2);
        end
        if y > size(filtCC1,1)
            y = size(filtCC1,1);
        end
        if z > size(filtCC1,3)
            z = size(filtCC1,3);
        end

        nucleiblue(y,x,z) = 255;

        nucloctable = [nucloctable; nucnum, x, y, z];
    end
    nucleiRGBstack(:,:,:,3) = nucleiblue;

    nucstackperm = permute(nucleiRGBstack, [1, 2, 4, 3]);
   %% Actin Analysis
   
    
    
    
    %% Actin loop
    for i=1:(size(S2vol,1))
        if tableS2(i,1) >= volumethreshmin2
            specind2 = cell2mat(S2idlist{i,1});
            x = specind2(:,1);
            y = specind2(:,2);
            z = specind2(:,3);
            for j = 1:(size(specind2,1))
                filtCC2(y(j),x(j),z(j)) = true;
            end
            actincount = actincount + 1;
            clear specind2;
        end
    end
       
    %% Labeling Found Nuclei
    totzval = size(nucstackperm, 4);
    zvalues = nucloctable(:,4);
    labeledstack = nucstackperm;
    for zsel = 1:totzval
        stacksel = nucstackperm(:,:,:,zsel);
        zindexes = find(zvalues == zsel);
        if ~isempty(zindexes)
            for i = 1:(size(zindexes,1))
                specind = zindexes(i);
                index = nucloctable(specind,1);
                xval = nucloctable(specind,2);
                yval = nucloctable(specind,3);
                stacksel = insertText(stacksel,[xval,yval], index, 'Fontsize', 10, 'TextColor', 'white', 'BoxOpacity', 0);
                stacksel = insertMarker(stacksel, [xval,yval], 'Color', 'white', 'Size', 3);
            end
            labeledstack(:,:,:,zsel) = stacksel;
            clear stacksel;
        end
    end
    
    %% Output folder

    chosdir = uigetdir;
    cd (chosdir);
    time = datetime('now');
    time = datestr(time);
    mkdir (time);
    cd (time);
    
    %% Exporting
    
    outputFileName_tiff = "functiontesting.tiff";
    outputTableName = "functiontesting_labeled.xls";
    
    for K=1:size(labeledstack,4)
       imwrite(labeledstack(:, :, :, K), outputFileName_tiff, 'WriteMode', 'append');
    end
    writetable(S1totalsort, outputTableName);
    
%    %%Extra Exporting
%             
%     raworiginal_name = 'raworiginal_image.tiff';
%     for K=1:size(stacktiff_thrcut,4)
%        imwrite(stacktiff_thrcut(:, :, :, K), raworiginal_name, 'WriteMode', 'append');
%     end
% 
%     rawthresh_name = 'rawthresh_image.tiff';
%     for K=1:size(threshstack,4)
%        imwrite(threshstack(:, :, :, K), rawthresh_name, 'WriteMode', 'append');
%     end
% 
%     thresh_nuclei = 'thresh_with_nucleifilt.tiff';
%     for K=1:size(nucstackperm,4)
%        imwrite(nucstackperm(:, :, :, K), thresh_nuclei, 'WriteMode', 'append');
%     end


 %end

