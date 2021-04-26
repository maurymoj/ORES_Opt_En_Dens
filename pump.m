function varargout = pump(varargin)
%pump Pump thermodynamic model. 
%   pump(h_in, h_s, eta_p) calculates h_out.
%
%   pump(h_in, h_s, h_out) calculates eta_p.
%=========================================================================%

if nargin == 6
   if strcmp(varargin{1},'h_in') && strcmp(varargin{3},'h_s') && strcmp(varargin{5},'eta_p')
      h_in = varargin{2};
      h_s = varargin{4};
      eta_p = varargin{6};
      h_out = h_in + (h_s - h_in)./eta_p;
      varargout{1} = h_out;
   elseif strcmp(varargin{1},'h_in') && strcmp(varargin{3},'h_s') && strcmp(varargin{5},'h_out') 
      h_in = varargin{2};
      h_s = varargin{4};
      h_out = varargin{6};
      eta_p = (h_s - h_in)./(h_out - h_in);
      varargout{1} = eta_p;
   end
end

end