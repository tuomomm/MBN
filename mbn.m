%MBN algorithm implementation for MATLAB.
%
%input:
%  N: Number of nodes
%  p: connection probability (in case of binomial PDF) or a vector of length N where p(i) is the probability of a node having i-1 inputs
%  weights: vector of length 16 determining the preference of each of the motifs
%  alpha: The importance of breaking of motifs in comparison to formation of motifs. Use alpha=1 for the performance as described in the article (formation and breaking of motifs are as important)
%  pow: The efficiency of the scoring. Use pow=inf for a performance as described in the article, 0 for a random network, and a value between 0 and inf for an intermediate efficiency
%  adapt: Is the adaptation in use. Use 1 for the performance as described in the article.
%
%output: 
%  M: Connectivity matrix, such that M_ij is the connection _from_ i _to_ j.
%
%Tuomo Maki-Marttunen, 2013-2016

function M = mbn(N,p,weights,alpha,pow,adapt)

if nargin < 6 || isempty(adapt)
    adapt = 0;
end

if nargin < 5 || isempty(pow)
    pow = 1;
end

if nargin < 4 || isempty(alpha)
    alpha = 1;
end

if nargin < 3 || isempty(weights)
    weights = zeros(1,16);
end

if length(weights) == 16
  weights = weights([1 2 4 6 5 9 3 8 14 7 12 10 13 11 15 16]); 
  %from now on, use the motifs are numbered (and labeled) as [(empty motif "-2"), (one-edge "-1"), (one bidirected edge "0"), (motifs "1", "2", ... , "13" as in numbered in Milo et al. 2002)]
end

if length(weights) == 13
    weights = [0 0 0 weights(:)'];
elseif length(weights) ~= 16
    disp('Give weights for motifs as a vector of length 13 or 16!');
end
weights = weights(:)';
if adapt
  adaptMtx = [0 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
              0 0 1 1 2 0 1 0 0 0 0 0 0 0 0 0;...
              0 0 0 0 0 2 0 0 0 2 0 0 0 0 0 0;...
              0 0 0 0 0 2 0 2 0 0 0 0 0 0 0 0;...
              0 0 0 0 0 1 0 1 0 1 0 1 0 0 0 0;...
              0 0 0 0 0 0 0 0 1 0 1 0 1 0 0 0;...
              0 0 0 0 0 0 0 2 0 2 0 0 0 0 0 0;...
              0 0 0 0 0 0 0 0 1 0 0 0 1 1 0 0;...
              0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0;...
              0 0 0 0 0 0 0 0 0 0 1 0 1 1 0 0;...
              0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0;...
              0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0;...
              0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0;...
              0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0;...
              0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;...
              0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
  weights = ((eye(16)-~~adaptMtx/N)\weights')';
end

if nargin < 2 || isempty(p)
    p = 0.2;
end

M = zeros(N);
if numel(p)==1  
  nns = binornd(N-1,p,N,1);  
elseif numel(p)==N
  nns = zeros(N,1);
  %p(i) gives the probability that number of in-neighbours is i-1 (i=1,...,N)
  p = p/sum(p);
  pcs = cumsum(p);
  for i=1:N
    [vain,nns(i)] = max(rand() <= pcs); %[~,nns(i)] not supported in all versions 
  end
  nns = nns - 1;
end

pinds = cell(N);
for i=1:N
  pinds{i} = [1:i-1,i+1:N];
end
nplaced = zeros(N,1);
  
%motifW: The number of motifs (-2 to 13) that are formed by drawing the
%connection from i through another node to k where first connection is
%bidi (mod 1), oppositely directed (mod 2), forward directed (mod 3) or
%non-connected (mod 0) and the connection from that to k is bidi (1-4),
%oppositedly (5-8), forward-directed (9-12) or nonexistent (13-16).
%On the condition that k _is not_ an output of i.
motifW = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;...%1
          0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0;...%9
          0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;...%5
          0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;...%13
          0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0;...%2
          0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;...%10
          0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;...%6
          0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0;...%14
          0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;...%3
          0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;...%11
          0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0;...%7
          0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;...%15
          0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;...%4
          0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0;...%12
          0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;...%8
          0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%16
%motifOW: The number of motifs (col -2 to 13) that are formed by drawing the
%connection from i through another node to k where first connection is
%bidi (mod 1), oppositely directed (mod 2), forward directed (mod 3) or
%non-connected (mod 0) and the connection from that to k is bidi (1-4),
%oppositedly (5-8), forward directed (9-12) or nonexistent (13-16).
%On the condition that k _is_ an output of i.
motifOW = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6;...%1
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;...%9
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;...%5
           0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0;...%13
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;...%2
           0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;...%10
           0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0;...%6
           0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;...%14
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;...%3
           0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0;...%11
           0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;...%7
           0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;...%15
           0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0;...%4
           0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;...%12
           0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;...%8
           0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0];%16
       
%motifminusW: The number of motifs (-2 to 13) that are broken by drawing the
%connection from i through another node to k where first connection is
%bidi (mod 1), oppositely directed (mod 2), forward directed (mod 3) or
%non-connected (mod 0) and the connection from that to k is bidi (1-4),
%oppositedly (5-8), forward directed (9-12) or nonexistent (13-16).
% On the condition that k _is not_ an output of i.
motifminusW = [0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0;...%1
               0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;...%9
               0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;...%5
               0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0;...%13
               0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;...%2
               0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;...%10
               0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0;...%6
               0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...%14
               0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;...%3
               0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0;...%11
               0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;...%7
               0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...%15
               0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0;...%4
               0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...%12
               0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...%8
               6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%16
%motifminusOW: The number of motifs (col -2 to 13) that are broken by drawing the
%connection from i through another node to k where first connection is
%bidi (mod 1), oppositely directed (mod 2), forward directed (mod 3) or
%non-connected (mod 0) and the connection from that to k is bidi (1-4),
%oppositedly (5-8), forward directed (9-12) or nonexistent (13-16).
% On the condition that k _is_ an output of i.
motifminusOW = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;...%1
                0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;...%9
                0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0;...%5
                0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;...%13
                0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;...%2
                0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0;...%10
                0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;...%6
                0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;...%14
                0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0;...%3
                0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;...%11
                0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;...%7
                0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0;...%15
                0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;...%4
                0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;...%12
                0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0;...%8
                0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%16
weightedmotifW = weights*( ~~motifW'-alpha*~~motifminusW');
weightedmotifOW = weights*(~~motifOW'-alpha*~~motifminusOW');
for ie=1:sum(nns) %go through all the edges in a non-obvious order

  w = nns - nplaced;
  w = w/sum(w);
  r = rand(1) < cumsum(w);
  i = find(r,1);
  
  pind = pinds{i};

  yesMandyesM = double(M&M');  %in and out
  yesMandnotM = double(M&~M'); %in but not out
  notMandyesM = double(~M&M'); %not in but out
  notMandnotM = double(~M&~M');%not in nor out
  notMandnotM = notMandnotM - eye(N); %The diagonal elements of M are zero, so the diagonal elements of notMandnotM would be ones.
                                      %These have to be cleared so that a node will not get extra non-neighbours.  
  
  %G is a 16-by-N matrix denoting the number of 2-paths from i to k where first is a bidi- (mod 1), in- (mod 2),
  %out- (mod 3) and non- (mod 0) connection and the second is a bidi (1-4), in- (5-8),
  %out- (9-12) and non- (13-16) connection 
  G = [reshape(yesMandyesM(i,:)*[yesMandyesM notMandyesM yesMandnotM notMandnotM],N,4)';...
       reshape(notMandyesM(i,:)*[yesMandyesM notMandyesM yesMandnotM notMandnotM],N,4)';...
       reshape(yesMandnotM(i,:)*[yesMandyesM notMandyesM yesMandnotM notMandnotM],N,4)';...
       reshape(notMandnotM(i,:)*[yesMandyesM notMandyesM yesMandnotM notMandnotM],N,4)'];
  %nNeighs is the number of bidi-, in-, out-, and non-neighbours of i
  %nNeighs = sum([yesMandyesM(:,i) yesMandnotM(:,i) notMandyesM(:,i) notMandnotM(:,i)]);
  points = M(i,pind).*(weightedmotifOW*G(:,pind)) + (~M(i,pind)).*(weightedmotifW*G(:,pind));

  if all(points <= 0) %if all points are negative, shift so that the minimum will be zero. this way
                      %the least negative will be dominant. if all were equally negative, the input
                      %by random. This step is needed only if pow is smaller than infinity
      points = points - min(points);
  else %otherwise set the negative weights to zero. since there are anyway some positive weights,
       %the input will be picked among those
      points(points<0) = 0;
  end
  if all(points==0)
     points = points+1;
  end
  if isinf(pow)
      points = points==max(points);
  else
      points = points.^pow;
  end
%  disp([num2str(points) ', inputs=' num2str(sum(M.*((1:N)'*ones(1,N))))]);
  points = points/sum(points);
  r = rand(1) < cumsum(points');
  rind = find(r,1);
  M(pind(rind),i) = 1;
  nplaced(i) = nplaced(i) + 1;
  pind(rind) = [];
  pinds{i} = pind;
end

