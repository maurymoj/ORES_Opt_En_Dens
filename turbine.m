function varargout = turbine(varargin)
%turbine    Turbine thermodynamic model.
%   turbine(h_in, h_s, eta_t)   calculates h_out
%   turbine(h_in, h_s, h_out)   calculates eta_t
%=========================================================================%

if nargin == 6
   if strcmp(varargin{1},'h_in') && strcmp(varargin{3},'h_s') && strcmp(varargin{5},'eta_t')
      h_in = varargin{2};
      h_s = varargin{4};
      eta_t = varargin{6};
      h_out = h_in - eta_t.*(h_in - h_s);
      varargout{1} = h_out;
   elseif strcmp(varargin{1},'h_in') && strcmp(varargin{3},'h_out') && strcmp(varargin{5},'h_s')
      h_in = varargin{2};
      h_out = varargin{4};
      h_s = varargin{6};
      eta_t = (h_in - h_out)./(h_in - h_s);      
      varargout{1} = eta_t;
   end
end

end