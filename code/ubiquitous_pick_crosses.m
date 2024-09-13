% ubiquitous project
% views ab-gamma crosses and user determines location

% Initialization
clc;
clear;
close all;

% Directory settings
figuresDir = 'E:\spectrolaminar\AttnData\core\cont\Crosses_taper_PreStim\'; % Directory where the figs are saved
figureFiles = dir(fullfile(figuresDir, '*.fig')); % Get all figure files

% Initialize variable to store user inputs paired with base file names
userInputs = struct('fileName', {}, 'inputs', {});

for i = 1:length(figureFiles)
    % Load and display the figure
    figFilePath = fullfile(figuresDir, figureFiles(i).name);
    open(figFilePath);
    
    % Pause to allow the user to view the figure
    disp(['Viewing figure: ', figureFiles(i).name]);
    pause(0.1); % Adjust the pause duration as needed
    
    % Prompt user for inputs
    prompt = 'Enter your inputs (e.g., [1, 2, 3]): ';
    userInput = input(prompt);
    
    % Store the inputs along with the base file name
    userInputs(end+1).fileName = figureFiles(i).name;
    userInputs(end).inputs = userInput;
    
    % Close the figure
    close(gcf);
end

% Save user inputs to a .mat file
save(fullfile(figuresDir, 'userInputsGammapks.mat'), 'userInputs');

disp('User inputs have been recorded and saved.');

% % Example user inputs data
% userInputsData = [10,18, 7,15, 11,16, 8,14, 3,12, 5,17, 3,17, 18, 2,15, ...
%     4,15, 8,14, 7,17, 6,13, 7,13, 8,14, 5,14, 6,11, 9,17, ...
%     7,17, 7,11, 12, 8,12, 4,14, 5,7,10,7,15, 10,17];
% 
% % Initialize an empty array to store the flattened inputs
% flattenedInputs = [];
% 
% % Loop through each entry in the userInputsData
% for i = 1:length(userInputs)
%     if iscell(userInputs{i})
%         % If the entry is a cell, convert it to a numeric array
%         flattenedInputs = [flattenedInputs, cell2mat(userInputs{i})];
%     else
%         % If the entry is a numeric array or single number, add it directly
%         flattenedInputs = [flattenedInputs, userInputs{i}];
%     end
% end
% 
% % Create a histogram of the user inputs
% figure;
% histogram(flattenedInputs, 'BinWidth', 1, 'FaceColor', 'b');
% title('Histogram of User Inputs');
% xlabel('Input Value');
% ylabel('Frequency');
% grid on;