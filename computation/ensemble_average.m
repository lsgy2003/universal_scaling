function avg_struct = ensemble_average(results)

    R = [results{:}];
    fields = fieldnames(R);
    avg_struct = struct();

    for i = 1:numel(fields)

        fname = fields{i};
        A = R(1).(fname);

        sz = size(A);
        tmp = zeros([sz, numel(R)]);

        for k = 1:numel(R)
            tmp(:,:,k) = R(k).(fname);
        end

        avg_struct.(fname) = mean(tmp, ndims(tmp));
    end
end