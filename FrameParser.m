classdef FrameParser

% Static methods
methods(Static)
    function frameStream = ParseFramesFromFiles(folder, ...
                                                generic_name, ...
                                                num_enumarator, ...
                                                second_name, ...
                                                extension)
        [~, s2] = size(num_enumarator);

        punc1 = '/';
        punc2 = '.';

        frameStream = [];
        
        for i = 1:s2
            r1 = imread(sprintf('%s%s%s%s%s%s%s', ...
                        folder, ...
                        punc1, ...
                        generic_name, ...
                        string(num_enumarator(i)), ...
                        second_name, ...
                        punc2, ...
                        extension));

            if size(frameStream,2) == 0
                frameStream = zeros(size(r1,1),size(r1,2),s2);
            end
                    
            r2 = r1(:,:,1:3);
            
            frameStream(:, :, i) = rgb2gray(r2);
        end
    end

    function frameStream = ParseFramesFromVideoReader(videoReader, ...
                                                      frames)
        % Initialize frame stream
        frameStream = zeros(videoReader.Height, videoReader.Width, size(frames,2));

        for i = 1:size(frames,2)
            % Set current time of the video reader to the desired frame
            videoReader.CurrentTime = i / videoReader.FrameRate;

            % Read current frame from the video reader
            frameStream(:,:,i) = rgb2gray(readFrame(videoReader));
        end
    end

    function intensityVector = GetIntensityVectorFromVideoReader(videoReader, ...
                                                                 point, ...
                                                                 frames)
        frame_count = size(frames,2);

        intensityVector = zeros(1,frame_count);

        for i = 1:frame_count
            current_frame = frames(i);
            % Set current time of the video reader to the desired frame
            videoReader.CurrentTime = current_frame / videoReader.FrameRate;

            current_frame = readFrame(videoReader);

            intensityVector(i) = current_frame(point(1), point(2));
        end
    end
    
    function [IntensityVectors, FileNames, XLimits] = ExportIntensityVectors(path)
        m_pathFiles = dir(path)

        FileNames = {};
        IntensityVectors = {};

        for i = 3:size(m_pathFiles,1)
            m_currentFilename = strcat(path,m_pathFiles(i).name);

            v = VideoReader(m_currentFilename);
            int_vector = FrameParser.GetIntensityVectorFromVideoReader(v,[100,100],1:v.FrameRate*v.Duration-1);
            
            FileNames{i-2} = filename;
            IntensityVectors{i-2} = int_vector;
        end
        
        XLimits = zeros(size(filenames,1)-2,2);
        for i = 1:size(IntensityVectors,2)
            plot(IntensityVectors{i})
            [x,y] = ginput(2);
            XLimits(i,1) = x(1);
            XLimits(i,2) = x(2);
        end

        save(strcat(path,'IntensityVectors'), ...
            'IntensityVectors', ...
            'FileNames', ...
            'XLimits');
    end
end

end