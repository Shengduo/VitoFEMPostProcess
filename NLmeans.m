function [output] = NLmeans(input,f,t,h,optMethod)
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 %  input: image to be filtered
 %  f: radius of similarity window
 %  t: radius of search window
 %  h: half-amplitude of noise
 %
 %  Original code authors: Jose Manjon Herrera & Antoni Buades (Mar 2006)
 %  Improved code author: Ivan Vlahinich & Vito Rubino, Caltech (Feb 2012)
 %      Added: average/regression option (Mar 2012)         
 %
 %  Implemented:
 %      Non-local 'NL-means' filter: 
 %          A. Buades, B. Coll and J.-M. Morel.
 %              "A non-local algorithm for image denoising"
 %      Non-local 'NL-means" filter improvement: 
 %          A. Buades, B. Coll, J.-M. Morel.
 %              "The staircasing effect in neighborhood filters"
 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % set default method if none given
 if nargin<=4
     optMethod = 'regression';  % choose: (default) 'regression' or 'average'
 elseif strcmpi(optMethod,'regression') || strcmpi(optMethod,'average')
 else error('Choose between two methods: "regression" or "average".');
 end
     
 % preliminary input parameters
 [m n]  = size(input);                          % size of input image
 output = zeros(m,n);                           % initialize output image
 input2 = padarray(input,[f f],'symmetric');    % pad edges of input
 kernel = make_kernel(f);                       % gauss kernel 
 kernel = kernel / sum(sum(kernel));            %   w/ some modification
 [ro,so]= ndgrid(1:2*f+1,1:2*f+1);              % base indices of
 indo   = sub2ind(size(input2),ro(:),so(:))-1;  %   a search window
 
 % start loop over all pixels
 for i=1:m
 for j=1:n
          
         % initialize indices and counts
         i1 = i+ f;                 % window center
         j1 = j+ f;                 %   in padded image

         rmin = max(i1-t,f+1);      % limits of
         rmax = min(i1+t,m+f);      %   image search area 
         smin = max(j1-t,f+1);
         smax = min(j1+t,n+f);

         lr  = rmax-rmin+1;         % number of matrices
         ls  = smax-smin+1;         %   to compare with 
         
    % assemble indices by adding to base search window
    rr  = repmat( (rmin-f:rmax-f) ,[1 ls]); 
    ss  = repmat( (smin-f:smax-f)',[1 lr]); ss = reshape(ss',1,[]);
    ind = sub2ind(size(input2), rr(:), ss(:));
    ind = repmat(indo,[1 ls*lr]) + repmat(ind',[(2*f+1)^2 1]);

    % compare all windows at once
    W1 = repmat(reshape(input2(i1-f:i1+f,j1-f:j1+f),[],1),[1 lr*ls]);
    W2 = input2(ind);
    w  = exp( -1/2/h^2 * (kernel(:)' * (W1-W2).^2) );
    z  = W2(ceil((2*f+1)^2/2),:);
    
    % set weight of source-source comparison to next highest
    ind_source    = sub2ind(size(input2),rr,ss) ... % find source index
                       == sub2ind(size(input2),i,j);
    w(ind_source) = -inf;                           % set weight to -inf
    w(ind_source) = max(w);                         % set weight to max

    % apply nL-meanas algorithm to update current pixel intensity
    switch lower(optMethod)
        case{'average'}
            sumw = sum(w);
            if sumw >= 0, output(i,j) = sum(w.*z) / sumw;
            else          output(i,j) = input(i,j); end
            
        case{'regression'}
            [x,y] = ind2sub(size(input2),ind(ceil((2*f+1)^2/2),:)); 
            abc   = lscov( [x' y' ones(size(x'))] , z' , w');
            output(i,j) = [x(ind_source) y(ind_source) 1] * abc;
    end
    
    % print loop progress
    loopProgress((i-1)*n+j,m*n,'remaining');
 end
 end
 
% note: not truly a gaussian kernel !
function [kernel] = make_kernel(f)              
 
kernel=zeros(2*f+1,2*f+1);   
for d=1:f    
  value= 1 / (2*d+1)^2 ;    
  for i=-d:d
  for j=-d:d  
    kernel(f+1-i,f+1-j)= kernel(f+1-i,f+1-j) + value ;
  end
  end
end
kernel = kernel ./ f;
 


        