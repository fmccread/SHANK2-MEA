function [idx1, idx2] = well_indicies(wells)

idx1 = [];
idx2 = [];

for i = 1:length(wells)

    well = wells(i);

    switch well

        case "A1"
            a = 1;
            b = 1;
        case "A2"
            a = 1;
            b = 2;
        case "A3"
            a = 1;
            b = 3;
        case "A4"
            a = 1;
            b = 4;
        case "B1"
            a = 2;
            b = 1;
        case "B2"
            a = 2;
            b = 2;
        case "B3"
            a = 2;
            b = 3;
        case "B4"
            a = 2;
            b = 4;
        case "C1"
            a = 3;
            b = 1;
        case "C2"
            a = 3;
            b = 2;
        case "C3"
            a = 3;
            b = 3;
        case "C4"
            a = 3;
            b = 4;
    end

idx1 = [idx1, a];
idx2 = [idx2, b];

end