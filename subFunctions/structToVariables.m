function [] = structToVariables(structIn)
try
	structFields 			= fieldnames(structIn);
	structIn_externalName 	= inputname(1);
	cellfun(@(fieldName) evalin('caller',[fieldName '=' structIn_externalName '.(' fieldName ');']), )
catch
	disp('Input is not a struct!');
end
end