classdef DHMUtility

% Static methods
methods(Static)

    function DepthsList = ReconstructPhaseSteppingVideos(IntensityVectors,FileNames,XLimits)
        DepthsList = {};

        for i = 1:size(IntensityVectors,2)
            Depths = DHMUtility.ReconstructPhaseSteppingVideo(IntensityVectors{i}, ...
                                                   string(FileNames{i}), ...
                                                   XLimits(i,:));
            DepthsList{i} = Depths;
        end
    end

    function Depths = ReconstructPhaseSteppingVideo(IntensityVector,FileName,XLimits)
        CurrentFileName = FileName
        CurrentVector = IntensityVector;

        v = VideoReader(CurrentFileName);
        last_frame = v.FrameRate * v.Duration;
        Frames = FrameParser.ParseFramesFromVideoReader(v,1:last_frame-1);

        y = CurrentVector(ceil(XLimits(1)):floor(XLimits(2)));
        x = 1:size(y,2);

        s = sinefit(x,y)
        period = floor(1/s(3));

        StartIndices = ceil(XLimits(1)):2:floor(XLimits(2))-period;
        Depths = zeros(v.Height,v.Width,size(StartIndices,2));

        counter = 1;
        for index = StartIndices
            CurrentDepth = Reconstructor.Reconstruct(double(Frames(:,:,index:index+period)),1);
            mesh(CurrentDepth);
            Depths(:,:,counter) = CurrentDepth;
            counter = counter + 1
        end
    end
    
    function ViewDepthVideo(Depths)
        for i = 1:size(Depths,3)
            mesh(Depths(:,:,i));
            pause(0.1);
        end
    end

    function PlotIntensityVector(video_reader,data_path,trimmed_name)
        % Extract frames
        frames = read(video_reader);
        frames = squeeze(frames(:,:,1,:));

        % Calculate the intensity vector
        point = [10 10];
        intensity_vector = squeeze(frames(point(1),point(2),:));
        
        frame_count = size(frames,3);
        
        % Fit a sinusoidal function over the intensity values
        x = 1:frame_count;
        s = sinefit(x,double(intensity_vector'));

        % Parameters of the sinusoidal function
        offset = s(1);
        amplitude = s(2);
        frequency = s(3);
        phase_shift = s(4);
        period = 1 / frequency;
        
        
        X = double(intensity_vector);
        Y = lowpass(X,20,1e3);
        y_plot = amplitude * sin(2*pi*frequency*x + phase_shift) + offset;

        figure('Position',[10 10 1920 1080]);
        plot(x,Y,'LineWidth',2); hold on; 
        plot(x,intensity_vector); hold on; 
        plot(x,y_plot,'LineWidth',2);
        legend(["Signal","LPF Signal","Sinusoidal"]);
        
        output_path = strcat(data_path,'Export\',trimmed_name,'_Plot.tiff');
        saveas(gcf,output_path);
    end
end

end