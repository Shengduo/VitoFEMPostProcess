function  out = loopProgress(jj,maxjj,arg)
% Purpose: Print loop progress to screen. 
% Variables:
%      jj  ... current loop iteration
%	 maxjj ... maximum loop iterations 		
%      arg ... 'elapsed'   - time / percent loop elapsed
%	           'remaining' - time / percent loop remaining
%      out ... string of what was printed on screen (e.g. to display on plot) 
% Note: When 'arg' not specified, 'remaining' is the default value
% Warning: Do not call twice inside a loop. Adds minimal overhead.
% Author: Ivan Vlahinich, Caltech

persistent lastCall lastOut startTime optArg

% are there enough input arguments?	
if(nargin>=2)
    
    % initialize variables on first iteration
    if(jj==1)
        lastCall  = int8(-1); 
        startTime = clock; 
        if(nargin<3), optArg='remaining'; else optArg=arg; end
    end

    % save time and skip some updates
    tmpCall = int8(jj/maxjj*100);
    if(lastCall~=tmpCall)

        % first step is unique	 	
        if(jj ~= 1)
            fprintf(1,repmat('\b',1,17));
        else
            fprintf(1,'%s',[optArg,': ']);
        end
      
        % elapsed/remaining time option
        lastCall    = tmpCall;	   
        elapsedTime = etime(clock,startTime);
        switch lower(optArg)
         case{'elapsed'}
          pp = lastCall;
          tt = elapsedTime;
         case{'remaining'}
          pp = int8(100) - lastCall;	
          tt = (maxjj-jj)/jj*elapsedTime;
         otherwise
          error('Error: -countLoop- can only report time elapsed or remaining.');
        end	
      
        % print to screen	
        tt = round(tt);
            HH = floor(tt/3600);
            MM = floor((tt-3600*HH)/60);
            SS = round(tt-60*MM-3600*HH);
        tt = sprintf('%02ih %02im %02is',HH,MM,SS); 
        fprintf(1,'%03i%%, %s', pp,tt);	
        
        % (optional)
        if nargout>0
            lastOut = [optArg, ': ' , sprintf('%03i%%, %s', pp,tt)];
        end
    end
    out = lastOut;
    
  % clean up after last step 
  if(jj==maxjj)
    fprintf(1,repmat('\b',1,19+length(optArg)));
    fprintf(1,'\b\n');
  end 

else
  error('Error: -countLoop- needs at least two input arguments...')
end
return