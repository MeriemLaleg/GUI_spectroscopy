function numVax = freadVAXG(fid, number, method)
% FREADVAXG Converts a number read in IEEE-LE format to its VAXG
% floating point representation.
%
% Usage:
%    fid = fopen('file.vaxg', 'r', 'ieee-le');
%    num = freadVAXG(fid, numElements, 'format')
%
% 'format' options are the same as FREAD.

switch method   
    case {'float32', 'single'}
      rawUINT32 = fread(fid,number,'uint32=>uint32');
      numVax = uint32le_to_VAXfloat( rawUINT32 );
      
    case {'float64', 'double'}
      rawUINT32 = fread(fid,2*number,'uint32=>uint32');%read 2 32bit numbers
      numVax = uint64le_to_VAXfloatG(rawUINT32);
      
    case {'float'}
      if intmax == 2147483647 %32bit OS float is 32 bits
	     rawUINT32 = fread(fid,  number,'uint32=>uint32');
         numVax = uint32le_to_VAXfloatG(rawUINT32);
      else
         rawUINT32 = fread(fid,2*number,'uint32=>uint32');%read 2 32bit numbers
         numVax = uint64le_to_VAXfloatG(rawUINT32);      
      end
      
    otherwise
      numVax    = fread(fid, number, method);

end

end