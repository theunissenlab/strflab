function out = getIDstring(in, outnum)

if ~exist('outnum', 'var')
    outnum = 30;
end

x = serialize(in);
x = int32(x);


x = [x 25:1107:10000];

x2 = int32(1:length(x))*2141;
x2 = mod(x2+sum(x), 6841) - mod(x2-length(x)+sum(diff(x(1:4:end)))*241, 3767);
x = mod(x+sum(x(2:2:end))*141221, 38021) - mod(x-sum(diff(x(1:3:end))), 758);
x = x - x2 + shift(x, 137) - shift(x, 17) + shift(x,3) - shift(x, 1);

sx = 0;
lx = length(x);
for ii=1:20
  tx = mod(sx, 10)+1;
  tx2 = mod(x(mod(sx,lx)+1),10)+1;
  sx = sx + sum(x(tx:tx2:end));
end
sx = mod(sx, lx)+1;

for ii=1:outnum
  sx = mod(sx + double(x(sx))+mod(ii*331, 521), lx)+1;
  out(ii) = mod(x(sx), 26)+65;
end

out = char(out);
