function handles=plot3vect(uuu,name,format,color,line_width)
%plot 3 vector  with tail at origin
%expects 
%   uuu=[ux uy uz]' (a column vector) 
%   name=string
%   format = plot format string (e.g., 'k:' for black dotted)
if nargin<4
    handles=plot3([0 uuu(1)]',[0 uuu(2)]',[0 uuu(3)]',format);
elseif nargin<5
    handles=plot3([0 uuu(1)]',[0 uuu(2)]',[0 uuu(3)]',format,'Color',color);
else
    handles=plot3([0 uuu(1)]',[0 uuu(2)]',[0 uuu(3)]',format,'Color',color,'LineWidth',line_width);
end
set(handles,'LineWidth',2)
text(uuu(1),uuu(2),uuu(3),name);