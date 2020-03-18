function [data]=open_binary(filename)

f=fopen(filename);
nfield=fread(f,1,'uint32');
data={};
for i=1:nfield
    data{i}.dims=fread(f,1,'uint32');
    data{i}.size=fread(f,data{i}.dims,'uint32')'; 
    data{i}.size=data{i}.size(end:-1:1); 
end
for i=1:nfield
    data{i}.data=reshape(fread(f,prod(data{i}.size),'double'),data{i}.size);
end
