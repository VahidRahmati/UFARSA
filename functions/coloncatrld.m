% This function is for "Vectorizing the Notion of Colon (:)", and was taken from (author: Loren Shure):
% http://blogs.mathworks.com/loren/2008/10/13/vectorizing-the-notion-of-colon/

function x = coloncatrld(start, stop)
% COLONCAT Concatenate colon expressions
%    X = COLONCAT(START,STOP) returns a vector containing the values
%    [START(1):STOP(1) START(2):STOP(2) START(END):STOP(END)].

% Based on Peter Acklam's code for run length decoding.
len = stop - start + 1;

% keep only sequences whose length is positive
pos = len > 0;
start = start(pos);
stop = stop(pos);
len = len(pos);
if isempty(len)
    x = [];
    return;
end

% expand out the colon expressions
endlocs = cumsum(len);  
incr = ones(1, endlocs(end));  
jumps = start(2:end) - stop(1:end-1);  
incr(endlocs(1:end-1)+1) = jumps;
incr(1) = start(1);
x = cumsum(incr);