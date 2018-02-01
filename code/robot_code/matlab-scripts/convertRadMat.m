function convertRadMat(inputFile, outputFile)
    S = load(inputFile);
    radTo1024 =@(x) round(x/2/pi*1024+512);
    A = structfun(radTo1024, S, 'UniformOutput', false);
    B = struct2cell(A);
    C = fieldnames(A);
    for i = 1:numel(B)
        dlmwrite(outputFile, C{i}, '-append');
        dlmwrite(outputFile, B{i}, '-append');
    end
end