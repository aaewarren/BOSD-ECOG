%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PROBABILISTIC LOCALISATION OF DEPTH ELECTRODE RECORDINGS

%Authored by Aaron Warren (aaron.warren@unimelb.edu)
	
%This MatLab script performs analysis steps described in the following paper:

%Intrinsic and secondary epileptogenicity in bottom-of-sulcus dysplasia
%Emma Macdonald-Laurs, Aaron E. L. Warren, Wei Shern Lee, Joseph Yuan-Mou Yang, Duncan MacGregor, Paul J. Lockhart, Richard J. Leventer, Andrew Neal, & A. Simon Harvey.
	
%Usage: run DepthElec_ProbLoc.m
	
%The following steps are performed:

%1. Load in a comma-separated table file containing patient/depth position IDs and the [x,y,z] coordinates of electrode "entry points" manually defined on brain pial surface (generated using Freesurfer). Note that there are multiple entry points/recordings per patient
%2. Load in BOSD sub-region drawings (.nii format), manually segmented from patients' T1-weighted MRI (i.e., bottom-of-sulcus [BOS] and top-of-sulcus [TOS])
%3. Determine all possible lines ("trajectories") between the depth electrode entry point and the coordinates of every voxel within the BOS sub-region
%4. Find coordinates of "candidate" depth electrode locations using pre-defined distances (in mm) between the entry point and target voxels from step 3 above. Pre-defined distances = distances between entry point and each channel of the depth electrode. Distances assume bipolar recordings (mid-point between adjacent channels)
%5. Assign each of the candidate depth electrode locations (from step 4) to one of the manually drawn BOSD sub-regions, by first finding the distance between the candidate location and the coordinates of every voxel within each sub-region (BOS or TOS), and then assigning to the sub-region where this distance is both: (i) smallest; and (ii) less than 1mm. If there are no distances that satisfy these two criteria, then the candidate position is not assigned to either of the sub-regions and instead is labelled as "outside".
%6. For each depth electrode channel, determine which of the sub-regions (including "outside" as a possibility) is most likely to have contained the recording, by labelling the channel by the name of the sub-region (or 'outside') to which the greatest number of candidate positions were assigned (in step 5). 
%7. Append the most probable channel location (BOS, TOS, or outside) to the table loaded in step 1 
	
%The following software dependencies are required: 
%Statistical parametric mapping (SPM) version 12: https://www.fil.ion.ucl.ac.uk/spm/software/spm12
%Resting-state fMRI Data Analysis Toolkit version 1.8: http://www.restfmri.net/forum/REST_V1.8

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
	
%add dependencies:
addpath /home/bin/spm12
addpath(genpath('/home/bin/REST_V1.8_130615'));

%load in table containing patient IDs and manually-defined [x,y,z] depth electrode entry coordinates (from Freesurfer):
T=readtable('/home/DepthElec_ProbLoc_Table.txt'); %Column 1 = patient ID/depth electrode entry ID; Columns 2-4 = x,y,z coordinates - see example table structure ("DepthElec_ProbLoc_Table_Example.txt")

%path to BOSD sub-region ROIs - these are binary masks drawn upon each patient's T1 image)
roidir='/home/BOSD_DRAWINGS'; %files are named *BOS.nii and *TOS.nii - prefaced with patient ID. For example, FCD001_BOS.nii

%define ROI lists, one including the label "outside"
roilist={'BOS','TOS'};
roilist_incoutside={'BOS','TOS','outside'};

%known depths of electrode contacts (from entry point) in mm - assuming these are bipolar recordings (halfway between contacts 1-2, 2-3, 3-4):
depths_bipolar=[13.6595 8.695 3.695]; %depths in mm from the entry point
	
%add some empty columns to the table - for saving most probable ("WIN") position of each recording
for bdepth=1:length(depths_bipolar);
    col1=append('BDEPTH_', num2str(depths_bipolar(bdepth)), 'MM_WIN');
    T.(col1)=cell(height(T),1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LOOP OVER PATIENTS/DEPTH ENTRIES DEFINED IN TABLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for d = 1:height(T);
 
    currdir=pwd;
    
    %extract patient ID and entry coordinates from table 
    subj=cellfun(@(s)s(1:6),T{d,1},'UniformOutput',false); %first six characters of first column = patient ID (e.g., FCD001). 
    pos=cellfun(@(s)s(7:10),T{d,1},'UniformOutput',false); %characters 7-10 = position ID for this patient (e.g., P08A). 
    entrycoordsRAS=T{d,2:4}; %x,y,z coords of depth electrode entry point (Freesurfer RAS format)
      
    disp(['starting analysis for ',char(subj),' depth position ',char(pos),' ...']);
    
    %load BOS ROI, find vox coords for non-zero voxels, and convert to mm
    targetROI=append(subj, '_BOS.nii'); %define target roi name, for importing
    Ytarg = spm_read_vols(spm_vol(char(fullfile(roidir, targetROI))),1);
    indx = find(Ytarg>0);
    [x,y,z] = ind2sub(size(Ytarg),indx);
    Tvoxcoords=[x,y,z];
    [Data, Header]=rest_ReadNiftiImage(char(fullfile(roidir, targetROI))); 
    Tmmcoords=zeros(size(Tvoxcoords));
    for q = 1:length(Tvoxcoords);
        x1=Tvoxcoords(q,1);
        y1=Tvoxcoords(q,2);
        z1=Tvoxcoords(q,3);
        temp_coords=Header.mat*[x1; y1; z1; 1];
        Tmmcoords(q,1)=temp_coords(1);
        Tmmcoords(q,2)=temp_coords(2);
        Tmmcoords(q,3)=temp_coords(3);
    end
    
    %calculate line/length/unit vector between each voxel in the target ROI & the entry point
    %then use this info to calculate coordinate for each depth electrode
    %contact with pre-defined distances along each of these target trajectories
    Tlines=zeros(size(Tmmcoords));
    Tlengths=zeros(size(Tmmcoords,1),1);
    Tvecs=zeros(size(Tmmcoords));
    depths_bipolar_coords=zeros(size(Tmmcoords,1),3,size(depths_bipolar,2)); %note: 3d matrix 
    for q = 1:length(Tmmcoords);
        Tlines(q,:)=Tmmcoords(q,:)-entrycoordsRAS; %line between target ROI voxel and entry point
        Tlengths(q)=sqrt(Tlines(q,1)^2 + Tlines(q,2)^2 + Tlines(q,3)^2); %length between target ROI voxel and entry point
        Tvecs(q,:)=Tlines(q,:)/Tlengths(q);
        parfor bdepth=1:length(depths_bipolar);
            %starting from entry coordinates, advance depth of contact
            %multiplied by unit vector to derive coordinates of candidate
            %contact position
            depths_bipolar_coords(q,:,bdepth)=entrycoordsRAS+depths_bipolar(bdepth)*Tvecs(q,:);
        end
    end
			    
    %for each candidate contact position calculated above, assign to one of the 2 possible ROIs: BOS or TOS
    %achieve this by first finding distance between the candidate contact position and every
    %voxel within each of the ROIs. 
    
    %create empty matrices to store the smallest distances (between each
    %pairing of ROI voxel and candidate contact position), at each depth, for each ROI 
    smallestdist2roivox_bipolar=zeros(length(depths_bipolar_coords),size(depths_bipolar,2),length(roilist));
    for roi=1:length(roilist);
        roiname=append(subj, '_', roilist{roi}, '.nii'); %define roi name 
        Y = spm_read_vols(spm_vol(char(fullfile(roidir, roiname))),1);
        indx = find(Y>0);
        [x,y,z] = ind2sub(size(Y),indx);
        vox_coords=[x,y,z];
        [Data Header]=rest_ReadNiftiImage(char(fullfile(roidir, roiname)));
        mm_coords=zeros(size(vox_coords));
        dist_bipolar=zeros(length(mm_coords), length(depths_bipolar_coords), length(depths_bipolar));
        for q = 1:length(vox_coords);
            x1=vox_coords(q,1);
            y1=vox_coords(q,2);
            z1=vox_coords(q,3);
            temp_coords=Header.mat*[x1; y1; z1; 1];
            mm_coords(q,1)=temp_coords(1);
            mm_coords(q,2)=temp_coords(2);
            mm_coords(q,3)=temp_coords(3);
                for coord=1:length(depths_bipolar_coords); 
                    for bdepth=1:length(depths_bipolar); 
                        dist_bipolar(q,coord,bdepth)=norm(mm_coords(q,:)-depths_bipolar_coords(coord, 1:3, bdepth));
                    end
                    for bdepth=1:length(depths_bipolar);
                        smallestdist2roivox_bipolar(coord,bdepth,roi)=min(dist_bipolar(:,coord,bdepth));
                    end 
                end      
        end
    end
        
    %now determine which roi each depth coord belongs to by finding the roi
    %containing the voxel with smallest distance away - recalling that this
    %minimum distance must be <1mm. if >1mm, then label as "outside". 
    depth_bipolar_roi_id=cell(size(Tmmcoords,1),size(depths_bipolar,2));
    for coord=1:size(Tmmcoords,1);
        for bdepth=1:size(depths_bipolar,2);
            mindists=smallestdist2roivox_bipolar(coord,bdepth,:);
            [M, I] = min(mindists);
            if M >=1; %%%%%if the minimum distance >=1mm ... 
                depth_bipolar_roi_id{coord,bdepth}='outside';
            else %%%if minimum distance is <1mm ... 
                depth_bipolar_roi_id{coord,bdepth}=roilist{I};
            end
        end
    end
    
    %lastly, for each depth, determine the roi (or "outside") with the greatest number of
    %candidate electrode positions assigned to it. 
	winroi_bipolar=cell(1,size(depths_bipolar,2));
    for bdepth=1:size(depths_bipolar,2);
        for roi=1:length(roilist_incoutside);
            idx=strfind(depth_bipolar_roi_id(:,bdepth), roilist_incoutside(roi));
            idx = find(not(cellfun('isempty', idx)));
            clear idx;
        end
        [s,~,j]=unique(depth_bipolar_roi_id(:,bdepth));
        winroi_bipolar{1,bdepth}=s{mode(j)};
    end
    
    %as final step, append most probable roi (or "outside") to the table
    T(d,5:7)=winroi_bipolar;
	    
    disp(['finished analysis for ',char(subj),' depth position ',char(pos),' hooray!']);
        
end
           