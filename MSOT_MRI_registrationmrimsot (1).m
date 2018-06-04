function MSOT_MRI_registrationmrimsot(MRI_file, MSOT_file, base_points, input_points)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Registration for MSOT/MRI
  % Image registration originally for MRI/FMT hybrid system for Vascularity project
  % Now this software is re-designed for registration of MSOT and MRI images
  % Wuwei Ren, Feb 2017
  % updated, 5,5,2017
  % version 2.0
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % MRI 
  % resolution: dx1,dy1,dz1
  % zeroZ: z_offset_1
  % MRI_raw
  % MRI_fixed
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % MSOT
  % resolution: dx2,dy2,dz2
  % zeroZ: z_offset_2
  % MSOT_raw
  % MSOT_trans
  
  %%

  addpath(genpath('../data/'));
  %% read MRI and MSOT datasets
  % read MRI
  if (isempty(MRI_file))
    [MRI_file_name,MRI_folder_name] = uigetfile({'*.tif';'*.tiff'},'Select MRI raw image stack in tif','C:\Users\renw\Desktop\data');
  else
    [MRI_file_name,MRI_folder_name] = imread(MRI_file);
  end
  info = imfinfo('../data/MRI image.tif');
  n_slice_MRI = numel(info);
  MRI_raw=zeros(256,256,n_slice_MRI);
  % mkdir(MRI_folder_name,'referenceMRI');
  for i = 1: n_slice_MRI;
      temp_img_rgb=imread('../data/MRI image.tif',n_slice_MRI+1-i);
      temp_img=temp_img_rgb(:,:,1);
      MRI_raw(:,:,i)=imresize(temp_img,[256 256]);
      tempt_output = zeros(256,256,3);
      tempt_layer_1=MRI_raw(:,:,i)./max(max(MRI_raw(:,:,i)));
      tempt_output(:,:,1) = tempt_layer_1;
      tempt_output(:,:,2) = tempt_layer_1;
      tempt_output(:,:,3) = tempt_layer_1;
      output_filename = ['reshaped_MRI_' num2str(i) '.tif'];
      imwrite(tempt_output,[MRI_folder_name 'referenceMRI\' output_filename]);
      figure,imagesc(MRI_raw(:,:,i));
      saveas(gcf, sprintf('../results/MRI_figure.png', i))
  end
  %%
  % read MSOT
  if (isempty(MSOT_file))
    [MSOT_file_name,MSOT_folder_name] = uigetfile({'*.tif';'*.tiff'},'Select MSOT raw image (structure) stack in tif',MRI_folder_name);
  else
    [MSOT_file_name,MSOT_folder_name] = imread(MSOT_file);
  end
  info = imfinfo(MSOT_file);
  n_slice_MSOT = numel(info);
  MSOT_raw=zeros(256,256,n_slice_MSOT);
  for i = 1: n_slice_MSOT;
      temp_img_rgb=imread(MSOT_file,i);
      temp_img=temp_img_rgb(:,:,1);
      MSOT_raw(:,:,i)=imresize(temp_img,[256 256]);
      MSOT_raw(240:256,:,i)=0;
      MSOT_raw(:,240:256,i)=0;
      figure,imagesc(MSOT_raw(:,:,i));
      saveas(gcf, sprintf('../results/MSOT_figure%02d.png', i))
  end
  
  %% parameter setting 
  % MRI 
  % resolution: 
  dx1 = 0.078125;%mm
  dy1 = 0.078125;%mm
  dz1=1.5;%mm
  thickness_1=1;%mm
  % zeroZ: 
  z_offset_1=0;%mm
  % x y z real coordinate
  x1 = dx1*(1:256);
  y1 = dy1*(1:256);
  z1 = dz1*(1:n_slice_MRI)+z_offset_1;
  % MRI_raw
  % MRI_fixed
  
  % MSOT
  % resolution: 
  dx2=0.1;%mm
  dy2=0.1;%mm
  dz2=0.3;%mm
  thickness_2=0.3;%mm
  % zeroZ: 
  z_offset_2=-3.3;%mm
  % x y z real coordinate
  x2 = dx2*(1:256);
  y2 = dy2*(1:256);
  z2 = dz2*(1:n_slice_MSOT)+z_offset_2;
  % MSOT_raw;
  % MSOT_trans;
  
  %% Step1: Interpolation and off-set
  % generating the same number of slices as MRI for MSOT
  MSOT_trans=zeros(256,256,n_slice_MRI);
  for i=1:n_slice_MRI
      containing_slice=0;
      for j=1:n_slice_MSOT
          distance_z1_z2=abs(z1(i)-z2(j));
          if distance_z1_z2<(thickness_1/2)
              containing_slice=containing_slice+1;
              MSOT_trans(:,:,i)=MSOT_trans(:,:,i) + MSOT_raw(:,:,j);
          end
      end
      if containing_slice>0
          MSOT_trans(:,:,i)=MSOT_trans(:,:,i)/containing_slice;
          disp(['slice #' num2str(i) ': find ' num2str(containing_slice) ' slices at position ' num2str(z1(i)) ' mm']);
      else 
          disp(['slice #' num2str(i) ': can not find MSOT images around MRI']);
      end
  end
  
  %% Step2: Find the transform matrix
  % select which slice you want to do transform
  % num_trans=input('select which slice you want to do transform:\n');
  num_trans=1
  if ((num_trans>0)&&(num_trans<=n_slice_MSOT))
      disp('transforming location identified\n');
  else
      disp('exceed!!!')
      return
  end
  
  float_image=MSOT_trans(:,:,num_trans);
  float_image=float_image./max(max(float_image));
  
  fixed_image=MRI_raw(:,:,num_trans);
  fixed_image=fixed_image./max(max(fixed_image));
  if (isempty(input_points) || isempty(base_points))
    h=cpselect(float_image,fixed_image);
    reply = input('base_points & input_points saved?\n','s');
    if (isempty(input_points))||(isempty(base_points))
      disp('wrong!!!!!!!!!!!!!!!!!!!!!!!!!!!')
      return
    end
  end
  
  % t = cp2tform(new,base_points,'nonreflective similarity');
  t = cp2tform(input_points,base_points,'nonreflective similarity');
  
  D = size(fixed_image);
  transed_MSOT_chosen = imtransform(float_image,t,'XData',[1 D(2)],'YData',[1 D(1)]);
  figure, imshowpair(fixed_image,transed_MSOT_chosen ,'montage')
  figure, imshowpair(fixed_image,transed_MSOT_chosen )
  
  
  %% Step3: Apply this transform to structural MSOT first
  transed_MSOT_output = zeros(256,256,n_slice_MRI);
  mkdir(MSOT_folder_name,'referenceMSOT');
  for i=1:n_slice_MRI
      
      float_image=MSOT_trans(:,:,i);
      float_image=float_image./max(max(float_image));
  
      transed_MSOT_output(:,:,i) = imtransform(float_image,t,'XData',[1 D(2)],'YData',[1 D(1)]);
      tempt_output = zeros(256,256,3);
      tempt_layer_1=transed_MSOT_output(:,:,i)./max(max(transed_MSOT_output(:,:,i)));
      tempt_output(:,:,1) = tempt_layer_1;
      tempt_output(:,:,2) = tempt_layer_1;
      tempt_output(:,:,3) = tempt_layer_1;
      output_filename = ['transformed_MSOT_' num2str(i) '.tif'];
      imwrite(tempt_output,[MSOT_folder_name 'referenceMSOT\' output_filename]);
  end
  % saveas(gcf, '../results/newfig.png')
  %% Step4: Apply this transform to other MSOT
  new_MSOT_option='N'
  % note: on local machine, uncomment following line to analyze further data
  % new_MSOT_option=input('Do you want to apply the same transform to a new MSOT dataset? Y/N [Y]: ', 's');
  while (new_MSOT_option=='Y')
      [MSOT_file_name,MSOT_folder_name] = uigetfile({'*.tif';'*.tiff'},'Select a new MSOT stack in tif',MRI_folder_name);
      info = imfinfo([MSOT_folder_name MSOT_file_name]);
      n_slice_MSOT = numel(info);
      MSOT_raw=zeros(256,256,n_slice_MSOT);
      for i = 1: n_slice_MSOT;
          temp_img_rgb=imread([MSOT_folder_name MSOT_file_name],i);
          temp_img=temp_img_rgb(:,:,1);
          MSOT_raw(:,:,i)=imresize(temp_img,[256 256]);
          MSOT_raw(240:256,:,i)=0;
          MSOT_raw(:,240:256,i)=0;
      end
      % MSOT_raw
      MSOT_trans=zeros(256,256,n_slice_MRI);
      for i=1:n_slice_MRI
          containing_slice=0;
          for j=1:n_slice_MSOT
              distance_z1_z2=abs(z1(i)-z2(j));
              if distance_z1_z2<(thickness_1/2)
                  containing_slice=containing_slice+1;
                  MSOT_trans(:,:,i)=MSOT_trans(:,:,i) + MSOT_raw(:,:,j);
              end
          end
          if containing_slice>0
              MSOT_trans(:,:,i)=MSOT_trans(:,:,i)/containing_slice;
              disp(['slice #' num2str(i) ': find ' num2str(containing_slice) ' slices at position ' num2str(z1(i)) ' mm']);
          else 
              disp(['slice #' num2str(i) ': can not find MSOT images around MRI']);
          end
      end
      transed_MSOT_output = zeros(256,256,n_slice_MRI);
      name_new_trans_MSOT =  input('Name the new transformed MSOT data: ', 's');
      mkdir(MSOT_folder_name,name_new_trans_MSOT);
      for i=1:n_slice_MRI
  
          float_image=MSOT_trans(:,:,i);
          float_image=float_image./max(max(float_image));
          transed_MSOT_output(:,:,i) = imtransform(float_image,t,'XData',[1 D(2)],'YData',[1 D(1)]);
          tempt_output = zeros(256,256,3);
          tempt_layer_1=transed_MSOT_output(:,:,i)./max(max(transed_MSOT_output(:,:,i)));
          tempt_output(:,:,1) = tempt_layer_1;
          tempt_output(:,:,2) = tempt_layer_1;
          tempt_output(:,:,3) = tempt_layer_1;
          output_filename = [name_new_trans_MSOT '_slice_#' num2str(i) '.tif'];
          imwrite(tempt_output,[MSOT_folder_name name_new_trans_MSOT '\' output_filename]);
      end    
      new_MSOT_option=input('Do you want to apply the same transform to a new MSOT dataset? Y/N [Y]: ', 's');
  end
  disp('Thanks, no more MSOT');
end