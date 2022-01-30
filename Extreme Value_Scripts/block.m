function c=block(data,blocksz,func),
%Calculates blocks' maxima or minima
%
%  USAGE: out=block(data,blocksz,func)
%
%   data: Data vector
%blocksz: Blocksize
%   func: Function to be applied to blocks, can be 'min' or 'max'
%
%    out: Blocks' maxima or minima
if strcmp(func,'max')|strcmp(func,'min')
nblocks=floor(length(data)/blocksz);
if nblocks*blocksz==length(data),
    datamat=reshape(data(1:nblocks*blocksz),blocksz,nblocks);
    exts=eval([func '(datamat,[],1)']);
else
        
    last=data((nblocks*blocksz)+1:end);
    datamat=reshape(data(1:nblocks*blocksz),blocksz,nblocks);
    exts=max(datamat,[],1);
    exts=eval([func '(datamat,[],1)']);
    exts=eval(['[exts'';' func '(last)]']);
    end
c=exts;    


else
    error('func input can be ''min'' or ''max''')

end

