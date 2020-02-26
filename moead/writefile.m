clear
clc

FID = fopen('MOEAD.out', 'wt');

for i = 1 : 10
    FUN = load(['FUN', int2str(i), '.out']);
    for j = 1 : 100
        fprintf(FID, '%f\t', FUN(j, 1));
        fprintf(FID, '%f', FUN(j, 2));
        fprintf(FID, '\n');
    end
end
fclose(FID);