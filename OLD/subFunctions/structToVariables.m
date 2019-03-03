function [] = structToVariables(structIn)
try
	structIn_externalName 	= inputname(1);
    filedNames              = fieldnames(structIn);
    for fieldID = 1:numel(filedNames)
        fieldName = filedNames{fieldID};
        evalin('caller',[fieldName '=' structIn_externalName '.(''' fieldName ''');']);
    end
catch
	disp('Input is not a struct!');
end
end