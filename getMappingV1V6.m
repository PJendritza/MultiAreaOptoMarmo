function [ mappOut ] = getMappingV1V6( mappNameIn ,type, sessionID)
%GETMAPPINGV1V6 get channel mapping for individual shanks in V1 and V6
%   Call the function like this: getMappingV1V6('V1_1B','raw',321) and it will return
%   a vector of channel numbers in the order of superficial->deep
%   if mappNameIn is 'all', the complete channel list with the corresponding
%   shank names will be the output

if sessionID < 386
    % mapping for V6
    mapping.V6_1A =  [20;22;3;9;18;28;32;30;15;1;13;11;26;24;7;5;2;16;4;6;23;25;10;12;29;27;14;8;31;21;17;19];
    mapping.V6_1B =  [113;127;115;117;104;106;121;123;110;108;125;119;112;102;98;100;99;101;116;122;97;107;111;109;128;114;126;124;105;103;120;118];
    mapping.V6_2B =  [49;63;51;53;40;42;57;59;46;44;61;55;48;38;34;36;35;37;52;58;33;43;47;45;64;50;62;60;41;39;56;54];
    mapping.V6_2A =  [84;86;67;73;82;92;96;94;79;65;77;75;90;88;71;69;66;80;68;70;87;89;74;76;93;91;78;72;95;85;81;83];
    % <- superficial --------------------------------------------------------------------------------------- deep -> %
    % mapping for V1
    mapping.V1_1B =  [49;63;51;53;40;42;57;59;46;44;61;55;48;38;34;36;35;37;52;58;33;43;47;45;64;50;62;60;41;39;56;54];
    mapping.V1_1A =  [20;22;3;9;18;28;32;30;15;1;13;11;26;24;7;5;2;16;4;6;23;25;10;12;29;27;14;8;31;21;17;19];
    % <- superficial --------------------------------------------------------------------------------------- deep -> %
    
elseif sessionID == 386
    % mapping for V6
    mapping.V6_1A =  [61;59;29;23;63;53;49;51;17;31;19;21;55;57;25;27;15;1;13;11;41;39;7;5;35;37;3;9;33;43;47;45];
    mapping.V6_1B =  [113;127;115;117;104;106;121;123;110;108;125;119;112;102;98;100;99;101;116;122;97;107;111;109;128;114;126;124;105;103;120;118];
    mapping.V6_2B =  [42;40;44;46;10;8;34;36;4;6;38;48;2;12;16;14;30;28;60;50;32;22;18;20;56;58;54;52;24;26;64;62];
    mapping.V6_2A =  [84;86;67;73;82;92;96;94;79;65;77;75;90;88;71;69;66;80;68;70;87;89;74;76;93;91;78;72;95;85;81;83];
    % <- superficial --------------------------------------------------------------------------------------- deep -> %
    % mapping for V1
    mapping.V1_1B =  [49;63;51;53;40;42;57;59;46;44;61;55;48;38;34;36;35;37;52;58;33;43;47;45;64;50;62;60;41;39;56;54];
    mapping.V1_1A =  [20;22;3;9;18;28;32;30;15;1;13;11;26;24;7;5;2;16;4;6;23;25;10;12;29;27;14;8;31;21;17;19];
    % <- superficial --------------------------------------------------------------------------------------- deep -> %
elseif sessionID > 386
    % mapping for V6
    mapping.V6_1A = [61;59;29;23;63;53;49;51;17;31;19;21;55;57;25;27;15;1;13;11;41;39;7;5;35;37;3;9;33;43;47;45];
    mapping.V6_1B = [106;104;108;110;74;72;98;100;68;70;102;112;66;76;80;78;94;92;124;114;96;86;82;84;120;122;118;116;88;90;128;126];
    mapping.V6_2B = [42;40;44;46;10;8;34;36;4;6;38;48;2;12;16;14;30;28;60;50;32;22;18;20;56;58;54;52;24;26;64;62];
    mapping.V6_2A = [125;123;93;87;127;117;113;115;81;95;83;85;119;121;89;91;79;65;77;75;105;103;71;69;99;101;67;73;97;107;111;109];
    % <- superficial --------------------------------------------------------------------------------------- deep -> %
    % mapping for V1
    mapping.V1_1B =  [42;40;44;46;10;8;34;36;4;6;38;48;2;12;16;14;30;28;60;50;32;22;18;20;56;58;54;52;24;26;64;62];
    mapping.V1_1A =  [61;59;29;23;63;53;49;51;17;31;19;21;55;57;25;27;15;1;13;11;41;39;7;5;35;37;3;9;33;43;47;45];
    
    % <- superficial --------------------------------------------------------------------------------------- deep -> %
else
    error('Incorrect sessionID for channel mapping')
end


if strcmp(mappNameIn, 'all') && strcmp(type, 'raw')
    shankNames = fieldnames(mapping);
    for ishank = 1:numel(shankNames)
        for ich = 1:numel(mapping.(shankNames{ishank}))
            
            if contains(shankNames{ishank}, 'V1')
                mappOut.ChnNr(1,mapping.(shankNames{ishank})(ich)+128) = mapping.(shankNames{ishank})(ich);
                mappOut.shankName{1,mapping.(shankNames{ishank})(ich)+128} = shankNames{ishank};
            else
                mappOut.ChnNr(1,mapping.(shankNames{ishank})(ich)) = mapping.(shankNames{ishank})(ich);
                mappOut.shankName{1,mapping.(shankNames{ishank})(ich)} = shankNames{ishank};
            end
            
        end
    end
else
    
    mappOut = mapping.(mappNameIn);
    
    if strcmp(type, 'raw') && contains(mappNameIn, 'V1')
        mappOut = mappOut+128; % channels for raw data are 1-128 V6 and 129-192 V1
    elseif strcmp(type, 'normal') || strcmp(type, 'raw')
    else
        error('Incorrect datatype for channel mapping')
    end
    
end

end

