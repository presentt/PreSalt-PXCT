function Nr = edensity(int16data, varargin)
    % calculate electron density (electrons/cubic angstrom) 
    % from integer data and metadata structure or scaling factors

    floatingdata = single(int16data);

    switch length(varargin)
        case 3
            high_cutoff = varargin{1};
            low_cutoff = varargin{2};
            factor_edensity = varargin{3};
        case 1
            metadataStruct = varargin{1};
            high_cutoff = metadataStruct.high_cutoff;
            low_cutoff = metadataStruct.low_cutoff;
            factor_edensity = metadataStruct.factor_edensity;
    end

    Nr = (floatingdata.*(high_cutoff - low_cutoff)/(2^16-1) + low_cutoff).*factor_edensity;

    Nr = uint16(round(Nr.*1e4)); % to make integers meaningful?
    % WARNING: cannot handle electron densities higher than 6.5535
    % electrons/cubic angstrom (typically nonphysical)
end