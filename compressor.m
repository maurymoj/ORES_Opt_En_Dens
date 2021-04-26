function varargout = compressor(varargin)
%compressor Compressor thermodynamic model. 
%   compressor(h_in, h_s, eta_c) calculates h_out.
%
%   compressor(h_in, h_s, h_out) calculates eta_c.
%=========================================================================%

if nargin == 6
   if strcmp(varargin{1},'h_in') && strcmp(varargin{3},'h_s') && strcmp(varargin{5},'eta_c')
      h_in = varargin{2};
      h_s = varargin{4};
      eta_c = varargin{6};
      h_out = h_in + (h_s - h_in)./eta_c;
      varargout{1} = h_out;
   elseif strcmp(varargin{1},'h_in') && strcmp(varargin{3},'h_s') && strcmp(varargin{5},'h_out') 
      h_in = varargin{2};
      h_s = varargin{4};
      h_out = varargin{6};
      eta_c = (h_s - h_in)./(h_out - h_in);
      varargout{1} = eta_c;
   end
end

end