function u=build_stack_table(metadata,files)
% Given a txt file listing directories of all downloaded stacks that should
% be included in the table, this function builds a table with the following
% variables for each stack:
% set - set number (one set = a collection of experiments)
% exp - experiment number (one exp = one slice)
% roi - the id-number of the region of interest in the slice (one roi = one stack)
% vessel origin - whether stack comes from tumor or non-tumor area
% file - directory of stack-file
% setname - same as set but string
% expname - same as exp but string
% roiname - same as roi but string
% HGCC - cell line
% perturbation - control or treatment
% dose - dose of treatment
% unit - unit of dose
% deltat - interval of image acquision
% frames - number of frames in the image stack
% t - total duration length of image stack
%
% Input parameters:
% metadata - metadata on all experiments
% files - list of stacks (for example, see filez.txt)
%
% Output parameters:
% u - stacktable
%
% @authors: Sven Nelander
% @date: 14082024

% If no list of stacks is provided, use the filez.txt to build table
if(nargin==1)
    t=readtable('filez.txt','ReadVariableNames',false,'delimiter','\t');
    t=t.Var1;
else
    t=readtable(files,'ReadVariableNames',false,'delimiter','\t');
    t=t.Var1;
end

% Extract set, exp and roi from the stack list by iterating through the
% stacks
data=[];
vessel_info = {};
for i=1:length(t)
    s=t{i};
    
    [i1,i2]=regexp(s,'Set.*exp');
    s2=s(i1:i2); 
    s2=regexprep(s2,'Set ','');
    s2=regexprep(s2,'\/exp','');
    my_set=str2num(s2);

    % If the stacks are for vessel quantification, add information on
    % whether the stack comes from a tumor or non-tumor area
    if(contains(s,'vasculature'))
        vessel_status = extractBetween(s,'Vessels_', '/exp');
        vessel_info(i) = vessel_status;
    end

    [i1,i2]=regexp(s,'exp.*stack');
    s2=s(i1:i2); 
    s2=regexprep(s2,'exp\_','');
    [i1,i2]=regexp(s2,'^\d*');
    my_exp=str2num(s2(i1:i2));
    
    [i1,i2]=regexp(s,'exp.*stack');
    s2=s(i1:i2); 
    s2=regexprep(s2,'^exp.*roi\_','');
    s2=regexprep(s2,'\_.*','');
    my_roi=str2num(s2);

    try
        my_set(1);
    catch
        my_set = metadata.set(my_exp);
    end
    
    data(i,:)=[my_set my_exp my_roi];
end

% Convert to a table
u=array2table(data,'VariableNames',{'set','exp','roi'});
u.file=t;

if(contains(s,'vasculature'))
u.vessel_origin = vessel_info';
end

u=sortrows(u);

setname={};
expname={};
roiname={};
key1={};
for i=1:height(u)
    setname{i,1}=['set' num2str(u.set(i))];
    expname{i,1}=['exp' num2str(u.exp(i))];
    roiname{i,1}=['roi' num2str(u.roi(i))];
    key1{i,1}=[setname{i,1} expname{i,1}];
end
u.setname=setname;
u.expname=expname;
u.roiname=roiname;


metasetname={};
key2={};
for i=1:height(metadata)
    metasetname{i,1}=['set' num2str(metadata.set(i))];
    key2{i,1}=[metasetname{i,1} metadata.experiment{i,1}];
end
metadata.setname=metasetname;


[~,lookup]=ismember(key1,key2);

u.HGCC=metadata.HGCC(lookup);

u.perturbation=metadata.perturbation(lookup);
u.dose=metadata.dose(lookup);
u.unit=metadata.unit(lookup);
u.delta_t=metadata.delta_t(lookup);
u.frames=metadata.frames(lookup);
u.t=metadata.t(lookup);

end



