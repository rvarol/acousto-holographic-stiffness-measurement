classdef Reconstructor

methods(Static)
    function out = Wavelength(data)
        persistent staticVar;

        if nargin
            staticVar = data;
        elseif size(staticVar,1) == 0
            staticVar = 5.21;
        end

        out = staticVar;
    end

    function out = NormalizedSpeedOfLight
        persistent staticVar;

        if nargin
            staticVar = data;
        elseif size(staticVar,1) == 0
            staticVar = 1.633;
        end

        out = staticVar;
    end

    function [depth, amplitude, phase] = Reconstruct(frames, cycles)
        [M, N, ~] = size(frames);
        
        % Resolve raw phase.
        [amplitude, raw_phase] = Reconstructor.Resolve(frames, cycles);

        % Unwrap phase.
        phase = Reconstructor.PhaseUnwrap(raw_phase);

        % Obtain thickness.
        thickness = Reconstructor.Phase2Depth_Transmittance( ...
                        phase, ...
                        Reconstructor.Wavelength, ...
                        Reconstructor.NormalizedSpeedOfLight);
        
        X = repmat([1:M]',1,N);
        rX = reshape(X',[M*N,1]);
        rY = repmat([1:N]',M,1);
        rZ = reshape(thickness',[M*N,1]);
        
        fitZ = fit([rX,rY],rZ,'poly22');
        calZ = feval(fitZ,[rX,rY]);
        zz = reshape(calZ,[N,M])';

        %mesh(thickness);
        %figure;
        %mesh(zz);

        thickness = thickness - zz;
        depth = thickness;
    end
    
    function unwrapped_phase = PhaseUnwrap(wrapped_phase)
        [M,N] = size(wrapped_phase);

        % Get poisson's value.
        poisson_val = Reconstructor.CalculatePoisson(wrapped_phase);

        % Apply discrete cosine transform.
        poisson_val_dct = dct2(poisson_val);
        % Create unwrapped phase coefficients. in discrete cosine domain.
        phase_dct = zeros(M, N);

        % Calculate unwrapped phase coefficients.
        for i = 1:M
            for j = 1:N
                if ((i ~= 1) || (j ~= 1))
                    phase_dct(i,j) = poisson_val_dct(i,j) / ...
                        (2 * (cos(pi * (i - 1) / M) + cos(pi * (j - 1) / N) - 2));
                else
                    phase_dct(i,j) = poisson_val_dct(i,j);
                end
            end
        end

        % Apply inverste discrete cosine transform.
        unwrapped_phase = idct2(phase_dct);
    end

    function P = CalculatePoisson(A)
        % Assign given matrix to poisson's for the sake of simplicity of
        % creating a matrix in the same size.
        P = A;

        [M, N] = size(A);

        DX = diff(A,1,2);
        DY = diff(A,1,1);
        WDXD = wrapToPi([zeros(1,size(A,2)); DY]);
        WDYD = wrapToPi([zeros(size(A,1),1) DX]);
        WDX = wrapToPi([DY; zeros(1,size(A,2))]);
        WDY = wrapToPi([DX zeros(size(A,1),1)]);
        P = (WDX - WDXD) + (WDY - WDYD);
    end
    
    function depth = Phase2Depth_Transmittance(phase, wavelength, normalized_speed_of_light)
        k_free = 2 * pi / wavelength;
        k_medium = normalized_speed_of_light * k_free;

        depth  = phase * (k_free * k_medium / (k_free - k_medium));
    end
    
    function depth = Phase2Depth_Reflectance(phase, wavelength)
        depth = ((phase / (2 * pi)) * wavelength) / 2;
    end
    
    function [amplitude, phase] = Resolve(frames, cycle_count)
        intensity_squared = (double(frames)).^2;

        % Get size of the frames matrix.
        [M, N, I] = size(frames);

        % |R|^2 + |O|^2. 
        dc_squared = sum(intensity_squared, 3) / I;

        dev_sub = repmat(dc_squared, [1 1 I]);

        % Calculate the deviations cos(alpha)*(|R|Ox) + sin(alpha)*(|R|*Oy)
        deviations = (intensity_squared - repmat(dc_squared, [1 1 I])) / 2;

        % Name |R|Ox as a, |R|Oy as b, cos(alpha) = Ka, sin(alpha) = Kb
        % then it can be converted to system of linear equations. And name
        % deviations as D.
        % Ka * a + Kb * b = D.
        % If all them are column vectors, Ka(i) * a + Kb(i) * b = D(i). Ka
        % denoting cos(alpha(i)), Kb sin(alpha(i)). alpha(i) is equal to 
        % 2 * pi * (i - 1) / I.

        % Then we apply Pseudo Inverse matrix technique, {which can be expressed 
        % as (A*A)^-1 A*b = X, A(coefficient matrix), b(result matrix), X(unknowns matrix)}. 
        % x1 and x2's are unknowns, A is Ka and Kb's, b is deviations.
        alpha = cycle_count * 2.0 * pi * (0:I-1) / I;
        A = [cos(alpha); sin(alpha)]';
        Aplus = (A' * A) \ A';

        a = sum((deviations .* repmat(permute(Aplus(1, :), [1 3 2]), [M N 1])), 3);
        b = sum((deviations .* repmat(permute(Aplus(2, :), [1 3 2]), [M N 1])), 3);

        % Since a = |R|Ox and b = |R|Oy, a^2 + b^2 = |R|^2 * |O|^2. And dc_squared
        % is |R|^2 + |O|^2. Then this creates second order equation system. By using these, 
        % object and reference wave amplitude can be calculated.

        % Calculate delta(which ought to be positive).
        delta = dc_squared.^2 - 4 * (a.^2 + b.^2);

        % Find solutions. Accept the one with lower intensity as object wave.
        amplitude = squeeze(sqrt((dc_squared - sqrt(delta)) / 2));
        phase = squeeze(atan2(b, a));
    end
end

end