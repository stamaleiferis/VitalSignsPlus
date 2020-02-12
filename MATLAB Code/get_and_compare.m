function[] = get_and_compare(input_file, input_ending)
hold off
[rr, tm_2] = ann2rr(input_file, input_ending);

[sig, Fs, tm] = rdsamp(input_file, 1);
[amp_test, i_test, d_test, b_test] = pan_tompkin(sig, Fs, 0);


plot(tm, sig, '-o', 'MarkerIndices', tm_2, 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red');
hold on
plot(tm, sig, 'b-s', 'MarkerIndices', i_test, 'MarkerFaceColor', 'green');
xlabel('Time in Seconds');
ylabel('ECG values');
title(input_file);
legend('MIT', 'Pan-Tompkin');
end
