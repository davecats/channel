function isfit=fftfit(x)

x=uint64(x);

while mod(x,2)==0
    x = bitsra(x,1);
end
isfit=(x==1 || x==3);