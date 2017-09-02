
function print_progress_string(current_value, total_number, name)

    percent = round(current_value/total_number * 100.0);

    print_every_percent = 0.05;
    increment = total_number * print_every_percent;
    
    if (mod(current_value, increment) == 0)

        progress = sprintf('%s: %i%% (%i / %i)', name, percent, current_value, total_number);
        disp(progress)
    end

end