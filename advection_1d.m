function [output,x] = advection_1d(filename,nam,outnum)

temp = hdf5info([filename '_' int2str(outnum)]);
name = [temp.GroupHierarchy.Groups(1).Name '/' nam];

output = hdf5read([filename '_' int2str(outnum)],name);

lower = temp.GroupHierarchy.Groups(1).Groups(1).Attributes(5).Value;
upper = temp.GroupHierarchy.Groups(1).Groups(1).Attributes(6).Value;
cells = temp.GroupHierarchy.Groups(1).Groups(1).Attributes(4).Value;

dx = (upper(1)-lower(1))/double(cells(1));
x = linspace(lower(1)+dx/2,upper(1)-dx/2,cells(1));
