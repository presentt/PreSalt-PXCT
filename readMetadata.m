function metadataStruct = readMetadata(filename)
    opts = delimitedTextImportOptions("NumVariables", 3);
    opts.DataLines = [1, 5];
    opts.Delimiter = ["#", "="];
    opts.VariableNames = ["Var1", "name", "value"];
    opts.SelectedVariableNames = ["name", "value"];
    opts.VariableTypes = ["string", "string", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts = setvaropts(opts, ["Var1", "name"], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["Var1", "name"], "EmptyFieldRule", "auto");
    raw_metadata = readtable(filename, opts);
    
    raw_metadata.Properties.RowNames = raw_metadata.name;

    metadataStruct.low_cutoff = raw_metadata.value("low_cutoff");
    metadataStruct.high_cutoff = raw_metadata.value("high_cutoff");
    metadataStruct.factor = raw_metadata.value("factor");
    metadataStruct.pixel_size = raw_metadata.value("pixel size");
    metadataStruct.factor_edensity = raw_metadata.value("factor_edensity");
end