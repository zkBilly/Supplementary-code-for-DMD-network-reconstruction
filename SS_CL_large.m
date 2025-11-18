clear all

T = 100000; %100000
time_step = 0.01;  % Time step for evolution
N = 100;
rng(20)

adjacency_matrix = load('.\SmallWorld_adjacency_matrix_k=4_p=0.2_lap.mat');
adjacency_matrix = adjacency_matrix.adjacency_matrix;
coupling_strength = 0.1;  % 0.1
Laplace_matrix = adjacency2laplace(adjacency_matrix);
Laplace_matrix = coupling_strength * Laplace_matrix;
[APL, L] = avgPathLength(adjacency_matrix);
[ACC, C] = avgClusteringCoefficient(adjacency_matrix);

Koopman_time_step = 0.35;  % 0.35
% If using previously generated data, comment out this section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set other parameters for CL model
a = 10+(10-10).*rand(N,1);
r = 20+(30-20).*rand(N,1);
b = 8/3+(8/3-8/3).*rand(N,1);
params = {N, a, r, b, Laplace_matrix};
CL_Trajectory = zeros(3*N, round(T/time_step));
CL_Trajectory(:, 1) = ones([3*N, 1]);
progress = 0;
for i = 2:round(T/time_step)
    if mod(i, round(T/time_step)/10) == 0
        progress = progress + 10;
        disp(['Simulation signal generation progress: ', num2str(progress), '%'])
    end
    CL_Trajectory(:, i) = Runge_Kutta(CL_Trajectory(:, i-1), params, time_step, @CL);
end
CL_Trajectory = CL_Trajectory';

delete_t = 100;
delete_step = delete_t / time_step;
CL_Trajectory = CL_Trajectory(delete_step+1:end, :);

figure(8)
plot_i = 75;
plot3(CL_Trajectory(1:end, plot_i), CL_Trajectory(1:end, plot_i+N), CL_Trajectory(1:end, plot_i+2*N))
hold on
scatter3(CL_Trajectory(1, plot_i), CL_Trajectory(1, plot_i+N), CL_Trajectory(1, plot_i+2*N), 'r*')
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Coupled Lorenz');
grid on; % Display grid

figure(9)
subplot(3,1,1)
plot(CL_Trajectory(end-1000:end, plot_i:plot_i))
title('x')
subplot(3,1,2)
plot(CL_Trajectory(end-1000:end, plot_i+N:plot_i+N))
title('y')
subplot(3,1,3)
plot(CL_Trajectory(end-1000:end, plot_i+2*N:plot_i+2*N))
title('z')

CL_Trajectory = CL_Trajectory(:, 1+2*N:3*N); % Only keep z-dimensional signal

% % Perform Fourier transform on Z-axis projection signal
% fs = 2 * pi / time_step;
% N = 10000;  % FFT length N, continue using this value
% n = 0 : N-1;
% spectrum_f = n * fs / N - fs / 2;
% spectrum = abs(fftshift(fft(CL_Trajectory(:, 1), N)));
% figure(10)
% stem(spectrum_f, spectrum, 'r')
% title('FFT Spectrum')
% xlabel('Angular Frequency')
% ylabel('Amplitude')
% xlim([-20, 20])

% Use Hilbert transform to calculate average frequency of all oscillators
% Before performing Hilbert transform, must remove DC component from signal
dc_offset = mean(CL_Trajectory); % Calculate DC component (mean)
CL_Trajectory = CL_Trajectory - dc_offset; % Subtract DC component from signal
delete_t = 20;
delete_step = delete_t / time_step;
hilbert_CL_Traj = hilbert(CL_Trajectory);
hilbert_CL_Traj = hilbert_CL_Traj(delete_step+1:end, :);
phase = angle(hilbert_CL_Traj(2:end, :) ./ hilbert_CL_Traj(1:end-1, :));
% Use Hilbert transform to process FHN membrane potential data, then only keep phase information and discard amplitude information
CL_Trajectory = exp(1i*angle(hilbert_CL_Traj));

figure(11)
plot(real(CL_Trajectory(1:1001, plot_i)))
hold on
% When performing Koopman analysis, increase the interval for selecting points
% (n+m)*Koopman_time_step should be greater than the longest oscillator period,
% 1/Koopman_time_step should be greater than 2*maximum frequency
jump_step = Koopman_time_step / time_step;
if jump_step > 1
    Koopman_indices = find(mod(1:size(CL_Trajectory, 1), jump_step) == 1);
    CL_Trajectory_Koopman = CL_Trajectory(Koopman_indices, :);
    CL_Trajectory = CL_Trajectory_Koopman;
    phase_Koopman = sumMatrixByStep(phase, jump_step);
    phase = phase_Koopman;
    clear CL_Trajectory_Koopman phase_Koopman
else
    Koopman_time_step = time_step;
end
figure(11)
plot(1:jump_step:1001, real(CL_Trajectory(1:ceil(1001/jump_step), plot_i)))
title('z')

save('.\CL_Trajectory.mat', 'CL_Trajectory')
save('.\phase.mat', 'phase')
clear hilbert_FHN_Traj
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

CL_Trajectory = load('.\CL_Trajectory.mat');
CL_Trajectory = CL_Trajectory.CL_Trajectory;
phase = load('.\phase.mat');
phase = phase.phase;

segment_num = 5;
shift_step = 1000; % 1000
n = 5;  % Observable dimension 5
m = 250000;  % time dimension 250000
CL_Trajectory = CL_Trajectory(1:(n+m+(segment_num-1)*shift_step), :);
phase = phase(1:(n+m+(segment_num-1)*shift_step-1), :);
for i = 1:segment_num
    Average_frequency = sum(phase(1+(i-1)*shift_step:n+m-1+(i-1)*shift_step, :)) / ((n+m-1) * Koopman_time_step);
    Average_frequency_all{i} = Average_frequency'; 
end
clear phase
% figure(2)
% histogram(Average_frequency, 0:2.5:50);
figure(2)
plot(real(CL_Trajectory(1:100, 13)))
hold on
plot(real(CL_Trajectory(1:100, 99)))
hold on
plot(real(CL_Trajectory(1:100, 68)))
hold on
plot(real(CL_Trajectory(1:100, 56)))
title('z')
legend('13', '99', '68', '56')
% figure(10)
% plot(real(CL_Trajectory(1:100, 1:10)))
% title('v')

%% Perform Koopman analysis on signal segment j and record
% 1) Phase difference phase_difference_choose_j_n
% 2) Amplitude difference amplitude_difference_choose_j_n
% 3) Selected average frequency Average_frequency_choose_j
for j = 1:segment_num
    j
    CL_seg = CL_Trajectory(1+(j-1)*shift_step:n+m+(j-1)*shift_step, :);
    Average_frequency = Average_frequency_all{j};

    % Take states of all oscillators to form (n*N)*m Hankel matrix, data length used is n+m
    H = zeros(n * N, m+1);
    for k = 1:m+1
        for p =1:N
            H(((p-1)*n+1):p*n, k) = CL_seg(k : (n + k - 1), p);
        end
    end
    clear CL_seg
    
    % Singular Value Decomposition SVD
    [U, S, V] = svd(H, 'econ');  % U is new basis, S*V' are coordinates
    clear H
    % Dimensionality reduction
    % truncation_singular_value = 1e-3;
    truncation_singular_value = 0;
    for truncation_order = 1:size(S, 1)
        if S(truncation_order, truncation_order) <= truncation_singular_value
            truncation_order = truncation_order - 1;
            break
        end
    end
    U = U(:, 1:truncation_order);
    S = S(1:truncation_order, 1:truncation_order);
    V = V(:, 1:truncation_order);
    
    % Calculate matrix representation of Koopman operator in basis U
    X = S * V(1:m, :)';
    Y = S * V(2:(m+1), :)';
    K = Y * pinv(X);  % XX' is close to singular matrix, so use pseudo-inverse
    clear S V
    
    % Selection method: Find eigenvalues and eigenvectors of Koopman operator
    [eigenfunction, eigenvalue] = eig(K);
    clear K
    eigenvalue = diag(eigenvalue);
    Frequency = angle(eigenvalue) / Koopman_time_step;  % Frequency obtained here is angular frequency
    % Find time-domain eigenfunctions of Koopman operator
    eigenfunction = U * eigenfunction;  % Same column is same frequency, same row is same time
    min_eigenvalue = 0.85;
    % figure(6)
    % stem(Frequency, abs(eigenvalue))
    % hold on
    % plot([min(Frequency), max(Frequency)], [min_eigenvalue, min_eigenvalue], 'r--')
    % % title('Mode frequency and eigenvalue magnitude')
    % xlabel('$\omega$', 'Interpreter', 'latex')
    % ylabel('$A$', 'Interpreter', 'latex')
    % title('$angular\ frequency\ and\ magnitude\ of\ eigenvalue$', 'Interpreter', 'latex')
    clear U
    % New: Since we consider stable non-decaying oscillation modes, need to eliminate eigenvalues too far from 1
    Frequency = Frequency(abs(eigenvalue) > min_eigenvalue);
    eigenfunction = eigenfunction(:, abs(eigenvalue) > min_eigenvalue);
    eigenvalue = eigenvalue(abs(eigenvalue) > min_eigenvalue);

    % Calculate phase difference and amplitude difference
    Average_frequency_choose_j = zeros(size(Average_frequency));
    for centernode_n = 1:N
        Average_frequency_n = Average_frequency(centernode_n, 1);
        
        % Amplitude selection method: Select mode with maximum amplitude from Koopman modes
        near_num = numel(Frequency);
        f_diffs = abs(Frequency - Average_frequency_n);
        [~, idx] = sort(f_diffs);
        closestIndices = idx(1:near_num);
        % closestFrequency = Frequency(closestIndices);
        amplitude_centernode_n = sum(abs(eigenfunction(((centernode_n-1)*n+1):centernode_n*n, :)), 1)./n;
        amplitude_centernode_n = amplitude_centernode_n(closestIndices);
        [~, max_amplitude_index] = max(amplitude_centernode_n);
        max_amplitude_index = closestIndices(max_amplitude_index);
        Average_frequency_eigenfunction_choose = eigenfunction(:, max_amplitude_index);
        Average_frequency_choose_j(centernode_n, 1) = Frequency(max_amplitude_index, 1);
        % figure(7)
        % stem(Frequency(closestIndices), amplitude_centernode_n)
        % text(Frequency(max_amplitude_index, 1), max(amplitude_centernode_n), ...
        %     sprintf('(%g, %g)', Frequency(max_amplitude_index, 1), max(amplitude_centernode_n)), ...
        %     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
        %     'FontSize', 10, 'Color', 'r');
        % xlabel('$\omega$', 'Interpreter', 'latex')
        % ylabel('$average(|\varphi_{\omega_i}^k|)$', 'Interpreter', 'latex')
        % if j == 1
        %     saveas(figure(7), ['./DMD_spectrum/',num2str(centernode_n),'.jpg']);
        % end

        % Observe the shape of eigenmode
        % figure(8)
        % plot(abs(Average_frequency_eigenfunction_choose))
    
        phase_difference_choose_j_n = zeros(1, N);
        amplitude_difference_choose_j_n = ones(1, N);
        for k =1 : N
            if k ~= centernode_n
                phase_difference_choose_j_n(1, k) = sum(angle(Average_frequency_eigenfunction_choose(((k-1)*n+1):k*n, 1) ./ ...
                    Average_frequency_eigenfunction_choose(((centernode_n-1)*n+1):centernode_n*n, 1))) / n;
                amplitude_difference_choose_j_n(1, k) = sum(abs(Average_frequency_eigenfunction_choose(((k-1)*n+1):k*n, 1) ./ ...
                    Average_frequency_eigenfunction_choose(((centernode_n-1)*n+1):centernode_n*n, 1))) / n;
            else
                continue
            end
        end
        amplitude_difference_choose_j_n = log10(amplitude_difference_choose_j_n);

        phase_difference_choose{j, centernode_n} = phase_difference_choose_j_n;
        amplitude_difference_choose{j, centernode_n} = amplitude_difference_choose_j_n;
    end
    Average_frequency_choose{j} = Average_frequency_choose_j;
end
save('.\Average_frequency.mat', 'Average_frequency')
save('.\Average_frequency_choose.mat', 'Average_frequency_choose')
save('.\amplitude_difference_choose.mat', 'amplitude_difference_choose')
save('.\phase_difference_choose.mat', 'phase_difference_choose')
save('.\Average_frequency_choose.mat', 'Average_frequency_choose')

%% Automated network structure discrimination algorithm
adjacency_matrix = load('.\SmallWorld_adjacency_matrix_k=4_p=0.2_lap.mat');
adjacency_matrix = adjacency_matrix.adjacency_matrix;
Average_frequency = load('.\Average_frequency.mat');
Average_frequency = Average_frequency.Average_frequency;
amplitude_difference_choose = load('.\amplitude_difference_choose.mat');
amplitude_difference_choose = amplitude_difference_choose.amplitude_difference_choose;
phase_difference_choose = load('.\phase_difference_choose.mat');
phase_difference_choose = phase_difference_choose.phase_difference_choose;

% Set minimum threshold for amplitude attenuation max_amplitude_difference, value range (-inf, 0), approaches 0 as coupling strength increases
mad = -0;  % If greater than this value, consider two oscillators synchronized
% Set maximum Euclidean distance max_distance for (dp, da) points showing mode propagation from three signal segments
mrms = 0.024;  % Maximum centroid distance variance, if greater than this value, consider no edge between two nodes
acd = 1.5;  % Amplitude attenuation clustering distance amplitude_cluster_distance, value range [0, inf]

inference_adjacency_matrix = zeros(N);
synchronize_pair = {};  % Detect synchronized pairs in network by examining amplitude attenuation
for centernode_n = 1:N
    % Stability screening: Plot (dp, da) obtained from three signal segments on same plane
    % figure(2)
    % scatter(phase_difference_choose{1, centernode_n}, amplitude_difference_choose{1, centernode_n}, 72, "red",'filled')
    % for n = 1:N
    %     phase_difference_choose_centernode_n = phase_difference_choose{1, centernode_n};
    %     amplitude_difference_choose_centernode_n = amplitude_difference_choose{1, centernode_n};
    %     text(phase_difference_choose_centernode_n(1, n), amplitude_difference_choose_centernode_n(1, n), ['  ', num2str(n)], 'Color', 'red')
    % end
    % hold on
    % scatter(phase_difference_choose{2, centernode_n}, amplitude_difference_choose{2, centernode_n}, 72, "blue",'o')
    % for n = 1:N
    %     phase_difference_choose_centernode_n = phase_difference_choose{2, centernode_n};
    %     amplitude_difference_choose_centernode_n = amplitude_difference_choose{2, centernode_n};
    %     text(phase_difference_choose_centernode_n(1, n), amplitude_difference_choose_centernode_n(1, n), ['  ', num2str(n)], 'Color', 'blue')
    % end
    % hold on
    % scatter(phase_difference_choose{3, centernode_n}, amplitude_difference_choose{3, centernode_n}, 72, "green",'o')
    % for n = 1:N
    %     phase_difference_choose_centernode_n = phase_difference_choose{3, centernode_n};
    %     amplitude_difference_choose
    %     text(phase_difference_choose_centernode_n(1, n), amplitude_difference_choose_centernode_n(1, n), ['  ', num2str(n)], 'Color', 'green')
    % end
    % legend('segment1', 'segment2', 'segment3')
    % title('Distribution of (da, dp) on plane')
    % xlabel('Phase difference dp')
    % set(gca, 'XTick', (-pi) : pi/2 : pi)
    % set(gca, 'XTicklabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
    % xlim([-pi, pi])
    % ylabel('Amplitude attenuation log(da)')
    % grid on
    % saveas(figure(2), ['.\result\', num2str(centernode_n), '.jpg'])
    % close(figure(2))
    % 
    % % New: Plot single point
    % plot_n = 1;
    % figure(8)
    % scatter(phase_difference_choose{1, centernode_n}(plot_n), amplitude_difference_choose{1, centernode_n}(plot_n), 72, "red",'filled')
    % text(phase_difference_choose{1, centernode_n}(plot_n), amplitude_difference_choose{1, centernode_n}(plot_n), ['  ', num2str(plot_n)], 'Color', 'red')
    % hold on
    % scatter(phase_difference_choose{2, centernode_n}(plot_n), amplitude_difference_choose{2, centernode_n}(plot_n), 72, "blue",'o')
    % text(phase_difference_choose{2, centernode_n}(plot_n), amplitude_difference_choose{2, centernode_n}(plot_n), ['  ', num2str(plot_n)], 'Color', 'blue')
    % hold on
    % scatter(phase_difference_choose{3, centernode_n}(plot_n), amplitude_difference_choose{3, centernode_n}(plot_n), 72, "green",'o')
    % text(phase_difference_choose{3, centernode_n}(plot_n), amplitude_difference_choose{3, centernode_n}(plot_n), ['  ', num2str(plot_n)], 'Color', 'green')
    % legend('segment1', 'segment2', 'segment3')
    % title('Distribution of (da, dp) on plane')
    % xlabel('Phase difference dp')
    % set(gca, 'XTick', (-pi) : pi/2 : pi)
    % set(gca, 'XTicklabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
    % xlim([-pi, pi])
    % ylim([-4, 0])
    % ylabel('Amplitude attenuation log(da)')
    % grid on

    % Synchronization test: If amplitude attenuation of mode propagation is too small, i.e., greater than mad, consider two nodes synchronized
    da_seg1_n = amplitude_difference_choose{1, centernode_n};
    for i = 1:N
        % If natural frequencies are too close causing synchronization (amplitude attenuation too small), add both to synchronized pair
        if amplitude_difference_choose{1, centernode_n}(i)>=mad ...
                && amplitude_difference_choose{2, centernode_n}(i)>=mad ...
                && amplitude_difference_choose{3, centernode_n}(i)>=mad ...
                && amplitude_difference_choose{4, centernode_n}(i)>=mad ...
                && amplitude_difference_choose{5, centernode_n}(i)>=mad
            % New: Must share same mode to be considered synchronized
            if Average_frequency_choose{1, 1}(centernode_n) == Average_frequency_choose{1, 1}(i) ...
                    && Average_frequency_choose{1, 2}(centernode_n) == Average_frequency_choose{1, 2}(i) ...
                    && Average_frequency_choose{1, 3}(centernode_n) == Average_frequency_choose{1, 3}(i) ...
                    && Average_frequency_choose{1, 4}(centernode_n) == Average_frequency_choose{1, 4}(i) ...
                    && Average_frequency_choose{1, 5}(centernode_n) == Average_frequency_choose{1, 5}(i)
                if centernode_n ~= i
                    synchronize_pair{end + 1} = [centernode_n, i];
                    da_seg1_n(i) = NaN;  % If two oscillators synchronized, no longer consider edge between them
                else
                    da_seg1_n(i) = NaN;  % If two oscillators synchronized, no longer consider edge between them
                end
            end
        end
    end

    % % Step 1 - Stability screening: Calculate maximum Euclidean distance max_distance_centroid from centroid for same index points on (dp, da) overlap plot for centernode_n using segment_num signal segments
    % % If less than md, consider average mode from center node received by branch node is stable, therefore consider edge exists between them
    % max_distance_centroid = zeros(1, N);
    % for i = 1:N
    %     points = [phase_difference_choose{1, centernode_n}(i), da_seg1_n(i)];
    %     for s = 2:segment_num
    %         points = [points; [phase_difference_choose{s, centernode_n}(i), amplitude_difference_choose{s, centernode_n}(i)]];
    %     end
    %    max_distance_centroid(i) = maxDistanceFromCentroid(points); 
    % end
    % 
    % for i = 1:N
    %     if i ~= centernode_n
    %         if max_distance_centroid(i) <= md
    %             inference_adjacency_matrix(centernode_n, i) = round(amplitude_difference_choose{1, centernode_n}(i), 1);
    %         else
    %             da_seg1_n(i) = NaN;
    %         end
    %     end
    % end
    
    % New: Step 1 - Stability screening: Calculate variance of distances from centroid for same index points on (dp,da) overlap plot for centernode_n using segment_num signal segments
    % If this variance is less than mcr, consider average mode from center node received by branch node is stable, therefore consider edge exists between them
    centroid_rms = zeros(1, N);
    for i = 1:N
        points = [phase_difference_choose{1, centernode_n}(i), da_seg1_n(i)];
        for s = 2:segment_num
            points = [points; [phase_difference_choose{s, centernode_n}(i), amplitude_difference_choose{s, centernode_n}(i)]];
        end
       centroid_rms(i) = rmsOfDistanceFromCentroid(points); 
    end
    
    % New: Save original indices corresponding to centroid_rms sorted in ascending order
    [centroid_rms_sort, centroid_rms_sort_indices] = sort(centroid_rms);
    centroid_rms_sort_result = [centroid_rms_sort; centroid_rms_sort_indices];

    for i = 1:N
        if i ~= centernode_n
            if centroid_rms(i) <= mrms
                inference_adjacency_matrix(centernode_n, i) = round(amplitude_difference_choose{1, centernode_n}(i), 1);
            else
                da_seg1_n(i) = NaN;
            end
        end
    end

    %↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓If only performing stability screening, comment this section↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
    % Step 2 - Amplitude screening: Further use clustering algorithm to screen nodes with smaller attenuation based on amplitude attenuation information
    current_cluster = [];
    current_cluster_size = 0;
    for j = 1:N
        former_size = length(current_cluster);
        for i = 1:N
            current_data = da_seg1_n(i);
            if isempty(current_cluster)
                current_cluster = [max(da_seg1_n)];
            else
                % Calculate distance between current data and current cluster center
                distance = abs(current_data - mean(current_cluster));
                % If distance <= clustering distance, add data to current cluster
                if distance <= acd && distance ~= 0
                    current_cluster = [current_cluster, current_data];
                end
            end
        end
        % When data volume is small, need additional scan to prevent omissions
        current_cluster = unique(current_cluster);
        current_cluster_size = length(current_cluster);
        if current_cluster_size == former_size
            break
        end
    end
    % Construct row centernode_n of adjacency matrix based on above two steps results
    inference_adjacency_matrix(centernode_n, :) = ismember(da_seg1_n, current_cluster);
    inference_adjacency_matrix(centernode_n, centernode_n) = 0;
    %↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑If only performing stability screening, comment this section↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑
end
inference_adjacency_matrix(isnan(inference_adjacency_matrix)) = 0;

% Step 3 - Symmetry screening: Since MS method tests influence of one node's motion on dynamic behavior of remaining nodes, result is unidirectional
% If network is undirected network, then make all unidirectional edges bidirectional in matrix
for i = 1:N
    for j = 1:N
        if inference_adjacency_matrix(i, j) ~= inference_adjacency_matrix(j, i)  % If inference matrix is asymmetric
            % Keep
            % inference_adjacency_matrix(i, j) = 1;
            % inference_adjacency_matrix(j, i) = 1;
            % Remove
            inference_adjacency_matrix(i, j) = 0;
            inference_adjacency_matrix(j, i) = 0;
        end
    end
end

% Step 4 - To meet needs of coupling strength inference algorithm, default synchronized nodes have edges between them
% Also, at node level, all edges of synchronized nodes should be shared
% for pair_i = synchronize_pair
%     pair = pair_i{1};
%     inference_adjacency_matrix(pair(1), pair(2)) = 1;
%     inference_adjacency_matrix(pair(2), pair(1)) = 1;
%     inference_adjacency_matrix(pair(1), :) = inference_adjacency_matrix(pair(1), :) + inference_adjacency_matrix(pair(2), :);
%     inference_adjacency_matrix(:, pair(1)) = inference_adjacency_matrix(:, pair(1)) + inference_adjacency_matrix(:, pair(2));
%     inference_adjacency_matrix(pair(1), :) = inference_adjacency_matrix(pair(1), :) ~= 0;
%     inference_adjacency_matrix(:, pair(1)) = inference_adjacency_matrix(:, pair(1)) ~= 0;
%     inference_adjacency_matrix(pair(2), :) = inference_adjacency_matrix(pair(1), :);
%     inference_adjacency_matrix(:, pair(2)) = inference_adjacency_matrix(:, pair(1));
%     inference_adjacency_matrix(pair(1), pair(1)) = 0;
%     inference_adjacency_matrix(pair(2), pair(2)) = 0;
% end

% Calculate inference accuracy
TP = 0;
FP = 0;
TN = 0;
FN = 0;
TP_node = zeros(1, N);
FP_node = zeros(1, N);
TN_node = zeros(1, N);
FN_node = zeros(1, N);
for i = 1:N
    for j = 1:N
        if i == j 
            continue
        else
            if adjacency_matrix(i, j) == 1 && inference_adjacency_matrix(i, j) == 1
                TP = TP + 1;
                TP_node(i) = TP_node(i) + 1;
            elseif adjacency_matrix(i, j) == 0 && inference_adjacency_matrix(i, j) == 1
                FP = FP + 1;
                FP_node(i) = FP_node(i) + 1;
            elseif adjacency_matrix(i, j) == 0 && inference_adjacency_matrix(i, j) == 0
                TN = TN + 1;
                TN_node(i) = TN_node(i) + 1;
            elseif adjacency_matrix(i, j) == 1 && inference_adjacency_matrix(i, j) == 0
                FN = FN + 1;
                FN_node(i) = FN_node(i) + 1;
            end
        end
    end
end
TP = TP / 2;
FP = FP / 2;
TN = TN / 2;
FN = FN / 2;
recall = TP / (TP + FN)
precision = TP / (TP + FP)
beta = 1;  % beta value for F-score, >1 means more emphasis on recall, i.e., hope for fewer missed target edges (FN)
F_beta_score = ((1+beta^2)*TP)/((1+beta^2)*TP+beta^2*FN+FP)
recall_node = TP_node ./ (TP_node + FN_node);
precision_node = TP_node ./ (TP_node + FP_node);
F_beta_score_node = ((1+beta^2)*TP_node)./((1+beta^2)*TP_node+beta^2*FN_node+FP_node);

% Plot inferred target network based on adjacency matrix
G1 = graph(adjacency_matrix);
deg1 = degree(G1);
edge_num = sum(deg1)/2;
% nSizes1 = deg1.^(4/5)+2;  % 'MarkerSize', nSizes1
G2 = graph(inference_adjacency_matrix);
deg2 = degree(G2);
% nSizes2 = deg2.^(4/5)+2;  % 'MarkerSize', nSizes2
common_color_limits = [min(min(deg1), min(deg2)), ...
                       max(max(deg1), max(deg2))];
nColors1 = deg1;
nColors2 = deg2;

figure(1);
subplot(1, 2, 1)
f1 = plot(G1, 'Layout', 'circle', 'NodeCData', nColors1); % 'force' layout can layout nodes based on connection forces between nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
node_index = 1:N;
label = arrayfun(@num2str, node_index, 'UniformOutput', false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colormap parula
clim(common_color_limits)
labelnode(f1, 1:N, label)
title('Target Network')
f1.NodeFontWeight = 'bold';
f1.NodeLabelColor = 'r';
colorbar
subplot(1, 2, 2)
f2 = plot(G2, 'Layout', 'circle', 'NodeCData', nColors2); % 'force' layout can layout nodes based on connection forces between nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_beta_score_label = arrayfun(@num2str, round(F_beta_score_node, 2), 'UniformOutput', false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colormap parula
clim(common_color_limits)
labelnode(f2, 1:N, F_beta_score_label)
title('Resulting Network')
f2.NodeFontWeight = 'bold';
f2.NodeLabelColor = 'r';
colorbar
% For directed graphs, comment the following part
% f2.ArrowSize = 15;
% f2.EdgeLabel = G2.Edges.Weight;
% f2.EdgeFontSize = 10;
% f2.EdgeLabelColor = 'black';

%%
%%%%%%%%%% Recalculate inference accuracy after merging synchronization groups %%%%%%%%%%
% New: Calculate PLV of signals within synchronization groups
synchronize_frequency_diff = [];
synchronize_PLV = [];
for i = 1:size(synchronize_pair, 2)
    synchronize_frequency_pair{i} = Average_frequency(synchronize_pair{i}, 1);
    synchronize_frequency_diff = [synchronize_frequency_diff, ...
        2*abs(synchronize_frequency_pair{i}(1)-synchronize_frequency_pair{i}(2))/ ...
        (synchronize_frequency_pair{i}(1)+synchronize_frequency_pair{i}(2))];
    synchronize_PLV = [synchronize_PLV, PLV(CL_Trajectory(:, synchronize_pair{i}(1)), ...
        CL_Trajectory(:, synchronize_pair{i}(2)))];  % New: Calculate PLV of synchronized pairs
end
Average_synchronize_frequency_diff = mean(synchronize_frequency_diff);
Average_synchronize_PLV = mean(synchronize_PLV);

% New: Calculate PLV between all signals
all_PLV = zeros(N);
for i = 1:N
    for j = 1:N
        if i == j
            continue
        else
            all_PLV(i, j) = PLV(CL_Trajectory(:, i), CL_Trajectory(:, j));
        end
    end
end
Average_all_PLV = mean(mean(all_PLV));

% New: Remove incorrectly judged synchronized pairs based on frequency
keep_index = [];
for sfp = synchronize_frequency_pair
    if abs(sfp{1}(1)-sfp{1}(2)) > 1 % Add extra insurance layer
        keep_index(end+1) = 0;
    else
        keep_index(end+1) = 1;
    end
end
keep_index = logical(keep_index);
synchronize_pair = synchronize_pair(keep_index);
synchronize_frequency_pair = synchronize_frequency_pair(keep_index);
synchronize_group = synPair2synGroup(synchronize_pair);
for i = 1:size(synchronize_group, 2)
    synchronize_frequency{i} = Average_frequency(synchronize_group{i}, 1);
end
merged_adjacency_matrix = mergeSynchronizedNodeGroups(adjacency_matrix, synchronize_group, true);
merged_inference_adjacency_matrix = mergeSynchronizedNodeGroups(inference_adjacency_matrix, synchronize_group, true);
mergelabel = generateMergeLabel(node_index, synchronize_group);
merge_N = size(merged_adjacency_matrix, 1);

% Calculate inference accuracy
TP_m = 0;
FP_m = 0;
TN_m = 0;
FN_m = 0;
TP_m_node = zeros(1, merge_N);
FP_m_node = zeros(1, merge_N);
TN_m_node = zeros(1, merge_N);
FN_m_node = zeros(1, merge_N);
for i = 1:merge_N
    for j = 1:merge_N
        if i == j 
            continue
        else
            if merged_adjacency_matrix(i, j) == 1 && merged_inference_adjacency_matrix(i, j) == 1
                TP_m = TP_m + 1;
                TP_m_node(i) = TP_m_node(i) + 1;
            elseif merged_adjacency_matrix(i, j) == 0 && merged_inference_adjacency_matrix(i, j) == 1
                FP_m = FP_m + 1;
                FP_m_node(i) = FP_m_node(i) + 1;
            elseif merged_adjacency_matrix(i, j) == 0 && merged_inference_adjacency_matrix(i, j) == 0
                TN_m = TN_m + 1;
                TN_m_node(i) = TN_m_node(i) + 1;
            elseif merged_adjacency_matrix(i, j) == 1 && merged_inference_adjacency_matrix(i, j) == 0
                FN_m = FN_m + 1;
                FN_m_node(i) = FN_m_node(i) + 1;
            end
        end
    end
end
TP_m = TP_m / 2;
FP_m = FP_m / 2;
TN_m = TN_m / 2;
FN_m = FN_m / 2;
recall_m = TP_m / (TP_m + FN_m)
precision_m = TP_m / (TP_m + FP_m)
beta = 1;  % beta value for F-score, >1 means more emphasis on recall, i.e., hope for fewer missed target edges (FN)
F_beta_score_m = ((1+beta^2)*TP_m)/((1+beta^2)*TP_m+beta^2*FN_m+FP_m)
recall_m_node = TP_m_node ./ (TP_m_node + FN_m_node);
precision_m_node = TP_m_node ./ (TP_m_node + FP_m_node);
F_beta_score_m_node = ((1+beta^2)*TP_m_node)./((1+beta^2)*TP_m_node+beta^2*FN_m_node+FP_m_node);

% Plot inferred target network based on synchronized adjacency matrix
G3 = graph(merged_adjacency_matrix);
deg3 = degree(G3);
edge_num_m = sum(deg3)/2;
% nSizes3 = deg3.^(4/5)+2;  % 'MarkerSize', nSizes3
G4 = graph(merged_inference_adjacency_matrix);
deg4 = degree(G4);
% nSizes4 = deg4.^(4/5)+2; % 'MarkerSize', nSizes4
common_color_limits = [min(min(deg3), min(deg4)), ...
                       max(max(deg3), max(deg4))];
nColors3 = deg3;
nColors4 = deg4;

figure(4);
subplot(1, 2, 1)
f3 = plot(G3, 'Layout', 'circle', 'NodeCData', nColors3); % 'force' layout can layout nodes based on connection forces between nodes
colormap parula
clim(common_color_limits)
labelnode(f3, 1:merge_N, mergelabel)
title('Target network after merging synchronized groups', Interpreter='latex')
f3.NodeFontWeight = 'bold';
f3.NodeLabelColor = 'r';
colorbar
% Set axis label font style
set(gca,'fontsize', 12);
subplot(1, 2, 2)
f4 = plot(G4, 'Layout', 'circle', 'NodeCData', nColors4); % 'force' layout can layout nodes based on connection forces between nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_beta_score_m_label = arrayfun(@num2str, round(F_beta_score_m_node, 2), 'UniformOutput', false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colormap parula
clim(common_color_limits)
labelnode(f4, 1:merge_N, F_beta_score_m_label)
title('Resulting network after merging synchronized groups', Interpreter='latex')
f4.NodeFontWeight = 'bold';
f4.NodeLabelColor = 'r';
colorbar
% Set axis label font style
set(gca,'fontsize', 12);

toc

function Rms = rmsOfDistanceFromCentroid(points)
    % points is a matrix containing n rows and 2 columns, each row represents coordinates (x, y) of a 2D point
    % Calculate centroid
    centroid = mean(points, 1); % Calculate centroid coordinates by column average
    % Calculate distance from each point to centroid
    distances = pdist2(points, centroid);
    % Find distance of farthest point from centroid
    Rms = rms(distances);
end

function synGroup = synPair2synGroup(synPair)
    % inputGroups: input cell array of arrays
    synGroup = {};

    while ~isempty(synPair)
        currentGroup = synPair{1};
        synPair(1) = [];

        % Find other groups intersecting with current group
        intersectingIndices = findIntersectingGroups(currentGroup, synPair);

        % Merge intersecting groups
        while ~isempty(intersectingIndices)
            currentGroup = unique([currentGroup, cell2mat(synPair(intersectingIndices))]);
            synPair(intersectingIndices) = [];
            intersectingIndices = findIntersectingGroups(currentGroup, synPair);
        end

        % Add current group to result
        synGroup{end + 1} = currentGroup;
    end
end

function mergeLabel = generateMergeLabel(node_index, nodeGroupsToMerge)
    % Generate node labels after merging synchronization groups
    % Merge node_index and nodeGroupsToMerge and remove elements in nodeGroupsToMerge that contain node_index elements

    % Initialize result array
    mergeLabel = cell(1, length(nodeGroupsToMerge) + length(node_index));

    % Process label elements
    for i = 1:length(node_index)
        if ~any(cellfun(@(x) ismember(node_index(i), x), nodeGroupsToMerge))
            mergeLabel{i} = num2str(node_index(i));
        end
    end

    % Process nodeGroupsToMerge elements
    for i = 1:length(nodeGroupsToMerge)
        str = ['{', strjoin(arrayfun(@num2str, nodeGroupsToMerge{i}, 'UniformOutput', false), ', '), '}'];
        mergeLabel{length(node_index) + i} = str;
    end

    % Remove empty elements
    mergeLabel = mergeLabel(~cellfun('isempty', mergeLabel));
end

function mergedAdjMatrix = mergeSynchronizedNodeGroups(originalAdjMatrix, nodeGroupsToMerge, uniOption)
    % originalAdjMatrix: original adjacency matrix
    % nodeGroupsToMerge: node groups to merge, each row is a node group
    % uniOption: logical value, true means merged adjacency matrix only indicates edge existence, elements only contain {0,1}
    % ..................false means coupling strength in merged adjacency matrix is sum of coupling strengths on corresponding edges in original matrix

    % Copy original adjacency matrix
    mergedAdjMatrix = originalAdjMatrix;

    % Iterate through each node group, merge nodes and preserve edge information
    for i = 1:size(nodeGroupsToMerge, 2)
        nodeGroup = nodeGroupsToMerge{i};

        % Merge nodes
        mergedNode_row = sum(mergedAdjMatrix(nodeGroup, :), 1);
        mergedNode_column = sum(mergedAdjMatrix(:, nodeGroup), 2);
        
        if uniOption == true
            % Set non-zero values to 1
            mergedNode_row = mergedNode_row ~= 0;
            mergedNode_column = mergedNode_column ~= 0;
        end

        % Add merged node to adjacency matrix
        mergedAdjMatrix = [mergedAdjMatrix, mergedNode_column];
        mergedAdjMatrix = [mergedAdjMatrix; [mergedNode_row, 0]];

        % Remove original nodes
        mergedAdjMatrix(nodeGroup, :) = [];
        mergedAdjMatrix(:, nodeGroup) = [];

        % Update node indices in unmerged synchronization groups
        for j = i+1:size(nodeGroupsToMerge, 2)
            nodeGroupsToMerge{j} = updateNodeIndices(nodeGroupsToMerge{j}, nodeGroup);
        end
    end
end

function intersectingIndices = findIntersectingGroups(group, groups)
    % Find indices of other groups intersecting with given group
    intersectingIndices = [];
    for i = 1:length(groups)
        if any(ismember(groups{i}, group))
            intersectingIndices = [intersectingIndices, i];
        end
    end
end

function summedMatrix = sumMatrixByStep(matrix, jump_step)
    % Determine dimensions of new matrix
    numRows = ceil(size(matrix, 1) / jump_step);
    newMatrix = zeros(numRows, size(matrix, 2));
    
    % Iterate through original matrix, sum by jump_step
    for i = 1:numRows
        startRow = (i - 1) * jump_step + 1;
        endRow = min(i * jump_step, size(matrix, 1));
        newMatrix(i,:) = sum(matrix(startRow:endRow,:), 1);
    end
    
    summedMatrix = newMatrix;
end

function updatedGroup = updateNodeIndices(group, removedNodes)
    % Update node indices in unmerged synchronization groups

    % Get mapping from old indices to new indices for nodes in unmerged synchronization groups after previous synchronization group merge
    map = zeros(1, length(group));
    for i = 1:length(group)
        map(i) = sum(removedNodes < group(i));
    end
    
    % Update node indices in unmerged synchronization groups
    updatedGroup = group - map;
end

function x_n = Runge_Kutta(x_n_1, params, dt, ode)
    k1 = ode(x_n_1, params);
    k2 = ode(x_n_1 + dt * 0.5 * k1, params);
    k3 = ode(x_n_1 + dt * 0.5 * k2, params);
    k4 = ode(x_n_1 + dt * k3, params);
    x_n = x_n_1 + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
end

function dxyz = CL(xyz0, params)
    % params contains {N, a, r, b, Laplace_matrix}
    % delta = [];
    N = params{1};
    a = params{2};
    r = params{3};
    b = params{4};
    Laplace_matrix = params{5};
    
    x0 = xyz0(1:N, 1);
    y0 = xyz0(N+1:2*N, 1);
    z0 = xyz0(2*N+1:3*N, 1);
    dx = zeros(N, 1);
    dy = zeros(N, 1);
    dz = zeros(N, 1);
    for i=1:N
        zj = Laplace_matrix(i, :) * z0;
        dx(i, 1) = a(i)*(y0(i)-x0(i));
        dy(i, 1) = r(i)*x0(i)-y0(i)-x0(i)*z0(i);
        dz(i, 1) = x0(i)*y0(i)-b(i)*z0(i)+zj;
    end
    dxyz = [dx; dy; dz];
end

function plv = PLV(phase_sig1, phase_sig2)
    % Calculate PLV value between sequences, sequences are e^{j*phase}
    Ntrials = max(size(phase_sig1));
    % compute PLV
    e = phase_sig1 ./ phase_sig2;
    plv = abs(sum(e)) / Ntrials;
end

