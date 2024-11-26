function plotEnergies(energies, times, methods, file_name)
    % Plots energy values for different methods over time and optionally saves as a PDF.
    %
    % Parameters:
    %   energies - Cell array where each element is a vector of energy values for a method.
    %   times - Cell array where each element is a vector of time values corresponding to the energy values.
    %   methods - String array with the names of the methods.
    %   save_as_pdf - (Optional) String, the filename to save the plot as a PDF.

    % Ensure the inputs are consistent
    assert(length(energies) == length(times), 'Number of energy arrays must match number of time arrays');
    assert(length(energies) == length(methods), 'Number of methods must match number of energy/time arrays');

    % Initialize the figure
    figure;
    hold on;

    % Plot each method
    for i = 1:length(methods)
        plot(times{i}, energies{i}, 'DisplayName', methods{i}, 'LineWidth', 1.5);
        hold on;
    end

    % Add labels, title, and legend
    xlabel('Time');
    ylabel('Energy');
    title('Energy Comparison Across Methods');
    legend('Location', 'best');
    grid on;

    % Release the hold
    hold off;

    % Save the figure as a PDF if a filename is provided
    if nargin == 4 && ~isempty(file_name)
        saveas(gcf, file_name, 'pdf'); % Save the figure as a PDF
        disp(['Plot saved as: ', file_name]);
    end
end
