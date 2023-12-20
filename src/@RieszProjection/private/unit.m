function [out, xunit] = unit(w,xunit)
%UNIT simple unit conversions

c0 = RieszProjection.c0;
switch xunit
    case 'nm'
        out = 2e9*pi*c0./w;
    case 'um'
        out = 2e6*pi*c0./w;
        xunit = '\mu m';
    case 'mm'
        out = 2e3*pi*c0./w;
    case 'm'
        out = 2*pi*co./w;
    case 'PHz'
        out = w*1e-15/(2*pi);
        xunit = 'PHz';
    case 'THz'
        out = w*1e-12/(2*pi);
        xunit = 'THz';
    case 'Hz'
        out = w/(2*pi);
        xunit = 'THz';
    otherwise
        out = w;
        xunit = 's^{-1}';
end
end