%function [nuccount, S1totalsort, labeledstack, RGB_stack_labeled] = nuclei_analysis_2(app, stacktiff_thrcut, threshstack, parameters1, connectivity)
    
%% starter

[stackorg, pathname, index] = uigetfile('.tif', 'Open CELL Image');
stacktiff = tiffreadVolume(stackorg);

stacktiff_thrcut = stacktiff(:,20:400,20:100, :);

parameters1 = [50,500,0.2,0.5]; %volumethreshmin1,volumethreshmax1, extentthreshmin1,extentthreshmax1
parameters2 = [50]; %volumethreshmin2

threshstack(:,:,:,1) = imbinarize(stacktiff_thrcut(:,:,:,1)); %nuclei
threshstack(:,:,:,2) = imbinarize(stacktiff_thrcut(:,:,:,2)); %actin

connectivity = 6;

%% Hole filling

    %Filling actin holes
    actin_stack = stacktiff_thrcut(:,:,:,2);
    
    %z values
    actin_stack(:,:,1) = imfill(actin_stack(:,:,1), 'holes');
    actin_stack(:,:,size(actin_stack,3)) = imfill(actin_stack(:,:,size(actin_stack,3)), 'holes');
    
    % y values
    actin_stack = permute(actin_stack, [3 2 1]);
    
    actin_stack(:,:,1) = imfill(actin_stack(:,:,1), 'holes');
    actin_stack(:,:,size(actin_stack,3)) = imfill(actin_stack(:,:,size(actin_stack,3)), 'holes');
    
    % x values
    actin_stack = permute(actin_stack, [1 3 2]);
    
    actin_stack(:,:,1) = imfill(actin_stack(:,:,1), 'holes');
    actin_stack(:,:,size(actin_stack,3)) = imfill(actin_stack(:,:,size(actin_stack,3)), 'holes');
    
    %reshape
    
    actin_stack = permute(actin_stack, [2 3 1]);
    stacktiff_thrcut(:,:,:,2) = actin_stack;
        


%% 

            %parameters for threshholding
            volumethreshmin1 = parameters1(1); 
            volumethreshmax1 = parameters1(2);
            extentthreshmin1 = parameters1(3);
            extentthreshmax1 = parameters1(4);
            
            % Hole Filling
        
            threshstack(:,:,:,1) = imfill(threshstack(:,:,:,1), 'holes');
            
            %Pre-analysis
        
            CC1 = bwconncomp(threshstack(:,:,:,1), connectivity);
        
        
        
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
            
            S1total = splitvars(S1total); %split variables
            
            S1sorted = sortrows(S1total,4); %sort table
            
            namelist1 = linspace(1,size(S1total,1), size(S1total,1));
            namelist1 = permute(namelist1, [2,1]);
            tablelist1 = table(namelist1);
            tablelist1.Properties.VariableNames(1) = "Number";
            S1totalsort = [tablelist1 S1sorted];
        
            nucleiRGBstack(:,:,:,1) = filtCC1 * 255; %binarized layer
            
            nucleiRGBstack(:,:,:,2) = im2double(stacktiff_thrcut(:,:,:,2));
            
            nucleiblue = zeros(size(filtCC1));
            
            
                %% Mosiac Initialization
    
            bin_filt_CC1 = bwconncomp(filtCC1);
            label_CC1 = labelmatrix(bin_filt_CC1);

            size_filt_x = size(label_CC1,1);
            size_filt_y = size(label_CC1,2);
            size_filt_z = size(label_CC1,3);

            RGB_labeled = zeros(size_filt_x, size_filt_y, 3, size_filt_z);

            for z_size = 1:size(label_CC1,3)
                RGB_labeled(:,:,:,z_size) = label2rgb(label_CC1(:,:,z_size),'jet', 'k');
            end
            %% Identifying nuclei, setting them to 255 in blue stack
            nucloctable = [];
            nuccount = size(S1totalsort,1);
            for i = 1:nuccount
                x = round(S1totalsort{i,3});
                y = round(S1totalsort{i,4});
                z = round(S1totalsort{i,5});
                nucnum = S1totalsort{i,1};
        
        
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
            
            %% Labeling Found Nuclei
            totzval = size(nucstackperm, 4);
            zvalues = nucloctable(:,4);
            labeledstack = nucstackperm;
            RGB_stack_labeled = RGB_labeled;
            
            for zsel = 1:totzval
                stacksel = nucstackperm(:,:,:,zsel);
                RGBsel = RGB_labeled(:,:,:,zsel);
                zindexes = find(zvalues == zsel);
                if ~isempty(zindexes)
                    for i = 1:(size(zindexes,1))
                        specind = zindexes(i);
                        index = nucloctable(specind,1);
                        xval = nucloctable(specind,2);
                        yval = nucloctable(specind,3);
                        stacksel = insertText(stacksel,[xval,yval], index, 'Fontsize', 10, 'TextColor', 'white', 'BoxOpacity', 0);
                        RGBsel = insertText(RGBsel,[xval,yval], index, 'Fontsize', 10, 'TextColor', 'white', 'BoxOpacity', 0);
                        stacksel = insertMarker(stacksel, [xval,yval], 'Color', 'white', 'Size', 3);
                    end
                    labeledstack(:,:,:,zsel) = stacksel;
                    RGB_stack_labeled(:,:,:,zsel) = RGBsel;
                    clear stacksel;
                    clear RGBsel;
                end
            end
        
%end