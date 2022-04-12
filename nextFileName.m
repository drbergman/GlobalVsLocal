function filename = nextFileName(folder,base_name,ndigs)

filename = sprintf('%s/%s%%0%dd.mat',folder,base_name,ndigs);
% filename = sprintf('data/%s%%03d.mat',base_name);
file_num = 1;
while exist(sprintf(filename,file_num),'file')
    file_num = file_num+1;
end

filename = sprintf(filename,file_num);