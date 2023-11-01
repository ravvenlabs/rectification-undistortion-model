clear all
load map.mat
for r = 1:480
  x = diff(map(r,:));
  map(r,1) = x(1);
  map(r,2:end) = x;
end

imshow(map)