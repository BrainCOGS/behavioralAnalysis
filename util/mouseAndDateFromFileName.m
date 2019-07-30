function [fmouse,fdate] = mouseAndDateFromFileName(fn)

% [fmouse,fdate] = mouseAndDateFromFileName(fn)
% deduces mouse name and date from file name

%% mouse name (try 4 characters first)
temp      = regexp(fn,'[_][a-z][a-z][0-9][0-9]','match');
if isempty(temp)
   temp  = regexp(fn,'[_][a-z][a-z][0-9][a-z]','match');
end
if isempty(temp)
   temp  = regexp(fn,'[_][a-z][a-z][0-9]','match');
end
if isempty(temp)
   temp  = regexp(fn,'[_][a-z][0-9][0-9]','match');
end
if isempty(temp)
   temp  = regexp(fn,'[_][a-z][0-9]','match');
end
if isempty(temp)
   temp  = regexp(fn,'[_][A-Z][0-9][0-9][0-9]','match');
end
if isempty(temp)
   temp  = regexp(fn,'[_][A-Z][0-9][0-9]','match');
end
if ~isempty(temp)
  fmouse = temp{1}(2:end);
else
  temp      = regexp(fn,'[a-z][a-z][0-9][0-9]','match');
  if isempty(temp)
    temp  = regexp(fn,'[a-z][a-z][0-9][a-z]','match');
  end
  if isempty(temp)
    temp  = regexp(fn,'[a-z][a-z][0-9]','match');
  end
  if isempty(temp)
    temp  = regexp(fn,'[a-z][0-9][0-9]','match');
  end
  if isempty(temp)
    temp  = regexp(fn,'[a-z][0-9]','match');
  end
  if isempty(temp)
    temp  = regexp(fn,'[A-Z][0-9][0-9][0-9]','match');
  end
  if isempty(temp)
    temp  = regexp(fn,'[A-Z][0-9][0-9]','match');
  end
  if isempty(temp)
    fmouse = [];
  else
    fmouse = temp{1};
  end
end

%% date
temp     = regexp(fn,'[0-9]{8,}','match');
if isempty(temp)
  fdate  = [];
else
  fdate  = temp{1};
end
