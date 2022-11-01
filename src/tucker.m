function T = tucker(T,varargin)
% TUCKER Tucker operator.
%    S = TUCKER(T, L) computes the Tucker operator
%
%       S = T x_1 L{1} x_2 L{2} x_3 ... x_d L{d}.
%
%    Here T is a complex tensor of size m_1 x ... x m_d, L a cell array
%    of complex matrices (L{mu} of size n_{mu} x m_{mu}) and x_mu denotes
%    the mu-mode product.
%
%    S = TUCKER(T, L1, L2, ..., Ld) computes the Tucker operator
%
%       S = T x_1 L1 x_2 L2 x_3 ... x_d Ld.
%
%    Here T is a complex tensor of size m_1 x ... x m_d, while Lmu is a
%    complex matrix of size n_{mu} x m_{mu}.
%
%    In both cases, if the entry corresponding to the mu-th matrix is empty,
%    then the associated mu-mode product is skipped.
%
%    [CCZ22] M. Caliari, F. Cassini, and F. Zivcovich,
%            A mu-mode BLAS approach for multidimensional tensor-structured
%            problems, Submitted 2022
  if (nargin < 2)
    error('Not enough input arguments.');
  end
  if (iscell (varargin{1}))
    varargin = varargin{1};
  end
  eidx = ~cellfun(@isempty, varargin);
  sT = [size(T), ones(1, length(varargin)-find(flip(eidx), 1)+1-length(size(T)))];
  lT = length(sT);
  mur = 1:length(varargin);
  mur = mur(eidx);
  lmu = length(mur);
  if (lmu == 0)
    error('Not enough non-empty input arguments.');
  end
  if (mur(1) == 1) && (mur(lmu) == lT)
    T = varargin{1}*reshape(T, sT(1), []);
    sT(1) = size(T, 1);
    T = reshape(T, [], sT(mur(lmu)))*varargin{mur(lmu)}.';
    sT(mur(lmu)) = size(T, 2);
    T = reshape(T, sT);
    mur = mur(2:lmu-1);
    lmu = lmu-2;
  elseif (mur(1) == 1) && (mur(lmu) ~= lT)
    T = varargin{1}*reshape(T, sT(1), []);
    sT(1) = size(T, 1);
    T = reshape(T, sT);
    mur = mur(2:lmu);
    lmu = lmu-1;
  elseif (mur(1) ~= 1) && (mur(lmu) == lT)
    T = reshape(T, [], sT(mur(lmu)))*varargin{mur(lmu)}.';
    sT(mur(lmu)) = size(T, 2);
    T = reshape(T, sT);
    mur = mur(1:lmu-1);
    lmu = lmu-1;
  end
  if (lmu > 0)
    T = permute(T, [mur(1), 1:(mur(1)-1), (mur(1)+1):lT]);
    for mu = 1:(lmu-1)
      T = varargin{mur(mu)}*reshape(T, sT(mur(mu)), []);
      sT(mur(mu)) = size(T, 1);
      T = permute(reshape(T, sT([mur(mu), 1:(mur(mu)-1), (mur(mu)+1):lT])), ...
      [mur(mu+1), 2:mur(mu), 1, (mur(mu)+1):(mur(mu+1)-1), (mur(mu+1)+1):lT]);
    end
    T = varargin{mur(lmu)}*reshape(T, sT(mur(lmu)), []);
    sT(mur(lmu)) = size(T, 1);
    T = ipermute(reshape(T, sT([mur(lmu), 1:(mur(lmu)-1), (mur(lmu)+1):lT])), ...
        [mur(lmu), 1:(mur(lmu)-1), (mur(lmu)+1):lT]);
  end
end