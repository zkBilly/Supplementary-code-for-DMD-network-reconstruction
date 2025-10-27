clear all

T = 100000; %100000
time_step = 0.01;  %演化的时间步长
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
%若使用先前生成的数据，注释此段%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 设置CL模型的其他参数
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
        disp(['仿真信号生成进度：', num2str(progress), '%'])
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
xlabel('X轴');
ylabel('Y轴');
zlabel('Z轴');
title('Couple Lorenz');
grid on; % 显示网格

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

CL_Trajectory = CL_Trajectory(:, 1+2*N:3*N); % 只保留z维信号

% %对Z轴投影信号作傅里叶变换
% fs = 2 * pi / time_step;
% N = 10000;  %FFT长度N，之后沿用此值
% n = 0 : N-1;
% spectrum_f = n * fs / N - fs / 2;
% spectrum = abs(fftshift(fft(CL_Trajectory(:, 1), N)));
% figure(10)
% stem(spectrum_f, spectrum, 'r')
% title('FFT频谱')
% xlabel('角频率')
% ylabel('幅度')
% xlim([-20, 20])

% 运用希尔伯特变换，计算所有振子的平均频率
% 进行希尔伯特变换之前，一定要去除信号中的直流分量
dc_offset = mean(CL_Trajectory); % 计算直流分量（均值）
CL_Trajectory = CL_Trajectory - dc_offset; % 从信号中减去直流分量
delete_t = 20;
delete_step = delete_t / time_step;
hilbert_CL_Traj = hilbert(CL_Trajectory);
hilbert_CL_Traj = hilbert_CL_Traj(delete_step+1:end, :);
phase = angle(hilbert_CL_Traj(2:end, :) ./ hilbert_CL_Traj(1:end-1, :));
% 运用希尔伯特变换处理FHN膜电位数据，进而只保留相位信息舍弃幅度信息
CL_Trajectory = exp(1i*angle(hilbert_CL_Traj));

figure(11)
plot(real(CL_Trajectory(1:1001, plot_i)))
hold on
%进行Koopman分析时，增大选点的间隔
% (n+m)*Koopman_time_step应大于最长的振子周期，
% 1/Koopman_time_step应大于2*最大频率
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
n = 5;  %可观测量维数指标（时间指标） 5
m = 250000;  %可观测量数量指标 250000
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

%% 对信号片段j进行Koopman分析，并记录：
% 1）相位差异phase_difference_choose_j_n
% 2）幅度差异amplitude_difference_choose_j_n
% 3）挑选出的平均频率Average_frequency_choose_j
for j = 1:segment_num
    j
    CL_seg = CL_Trajectory(1+(j-1)*shift_step:n+m+(j-1)*shift_step, :);
    Average_frequency = Average_frequency_all{j};

    %取所有振子的状态构成(n*N)*m的Hankel矩阵，所用的数据长度为n+m
    H = zeros(n * N, m+1);
    for k = 1:m+1
        for p =1:N
            H(((p-1)*n+1):p*n, k) = CL_seg(k : (n + k - 1), p);
        end
    end
    clear CL_seg
    
    %奇异值分解SVD
    [U, S, V] = svd(H, 'econ');  %U是新的基，S*V'是坐标
    clear H
    %降维
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
    
    %计算Koopman算符在基矢U下的矩阵表示
    X = S * V(1:m, :)';
    Y = S * V(2:(m+1), :)';
    K = Y * pinv(X);  %XX'接近奇异矩阵，因此用伪逆
    clear S V
    
    %挑选法：求Koopman算符的本征值和本征向量
    [eigenfunction, eigenvalue] = eig(K);
    clear K
    eigenvalue = diag(eigenvalue);
    Frequency = angle(eigenvalue) / Koopman_time_step;  %这里求得的频率是角频率
    %求Koopman算符的时域本征函数
    eigenfunction = U * eigenfunction;  %同一列是同一频率，同一行是同一时刻
    min_eigenvalue = 0.85;
    % figure(6)
    % stem(Frequency, abs(eigenvalue))
    % hold on
    % plot([min(Frequency), max(Frequency)], [min_eigenvalue, min_eigenvalue], 'r--')
    % % title('模式频率及本征值幅值')
    % xlabel('$\omega$', 'Interpreter', 'latex')
    % ylabel('$A$', 'Interpreter', 'latex')
    % title('$angular\ frequency\ and\ magnitude\ of\ eigenvalue$', 'Interpreter', 'latex')
    clear U
    %新增：因为我们考虑的是稳定不衰减的振荡模式，因此需要剔除本征值中过于远离1的值
    Frequency = Frequency(abs(eigenvalue) > min_eigenvalue);
    eigenfunction = eigenfunction(:, abs(eigenvalue) > min_eigenvalue);
    eigenvalue = eigenvalue(abs(eigenvalue) > min_eigenvalue);

    % 计算相位差异和幅度差异
    Average_frequency_choose_j = zeros(size(Average_frequency));
    for centernode_n = 1:N
        Average_frequency_n = Average_frequency(centernode_n, 1);
        
        % 幅度挑选法：从Koopman模式中挑选幅度最大的模式
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

        % 观察本征模式的形态
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

%% 自动化判别网络结构算法
adjacency_matrix = load('.\SmallWorld_adjacency_matrix_k=4_p=0.2_lap.mat');
adjacency_matrix = adjacency_matrix.adjacency_matrix;
Average_frequency = load('.\Average_frequency.mat');
Average_frequency = Average_frequency.Average_frequency;
amplitude_difference_choose = load('.\amplitude_difference_choose.mat');
amplitude_difference_choose = amplitude_difference_choose.amplitude_difference_choose;
phase_difference_choose = load('.\phase_difference_choose.mat');
phase_difference_choose = phase_difference_choose.phase_difference_choose;

% 设定幅度衰减的最小阈值max_amplitude_difference，取值范围为(-inf, 0)，依据耦合强度增大而趋近于0
mad = -0;  %大于此值则认为两振子发生了同步
% 设定从三段信号求得的展现某一模式传播情况的(dp, da)点的最大欧氏距离max_distance
mrms = 0.024;  %最大质心距离方差，大于此值则认为两节点间不存在连边
acd = 1.5;  %幅度衰减的聚类距离amplitude_cluster_distance，取值为[0, inf]

inference_adjacency_matrix = zeros(N);
synchronize_pair = {};  %通过检验振幅的衰减情况来探明网络中的同步对
for centernode_n = 1:N
    % 稳定性筛选：将三段信号得到的(dp, da)图画在同一平面上
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
    %     amplitude_difference_choose_centernode_n = amplitude_difference_choose{3, centernode_n};
    %     text(phase_difference_choose_centernode_n(1, n), amplitude_difference_choose_centernode_n(1, n), ['  ', num2str(n)], 'Color', 'green')
    % end
    % legend('segment1', 'segment2', 'segment3')
    % title('(da, dp)在平面上的分布')
    % xlabel('相位差dp')
    % set(gca, 'XTick', (-pi) : pi/2 : pi)
    % set(gca, 'XTicklabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
    % xlim([-pi, pi])
    % ylabel('振幅衰减log(da)')
    % grid on
    % saveas(figure(2), ['.\result\', num2str(centernode_n), '.jpg'])
    % close(figure(2))
    % 
    % % 新增：绘制单点
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
    % title('(da, dp)在平面上的分布')
    % xlabel('相位差dp')
    % set(gca, 'XTick', (-pi) : pi/2 : pi)
    % set(gca, 'XTicklabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
    % xlim([-pi, pi])
    % ylim([-4, 0])
    % ylabel('振幅衰减log(da)')
    % grid on

    % 同步检验：如果模式传播的幅度衰减程度过小，表现为大于mad，则认为两节点同步
    da_seg1_n = amplitude_difference_choose{1, centernode_n};
    for i = 1:N
        % 如果固有频率过于接近导致同步发生（幅度衰减过小），则将二者加入同步对
        if amplitude_difference_choose{1, centernode_n}(i)>=mad ...
                && amplitude_difference_choose{2, centernode_n}(i)>=mad ...
                && amplitude_difference_choose{3, centernode_n}(i)>=mad ...
                && amplitude_difference_choose{4, centernode_n}(i)>=mad ...
                && amplitude_difference_choose{5, centernode_n}(i)>=mad
            % 新增：必须共用同一模式才算同步
            if Average_frequency_choose{1, 1}(centernode_n) == Average_frequency_choose{1, 1}(i) ...
                    && Average_frequency_choose{1, 2}(centernode_n) == Average_frequency_choose{1, 2}(i) ...
                    && Average_frequency_choose{1, 3}(centernode_n) == Average_frequency_choose{1, 3}(i) ...
                    && Average_frequency_choose{1, 4}(centernode_n) == Average_frequency_choose{1, 4}(i) ...
                    && Average_frequency_choose{1, 5}(centernode_n) == Average_frequency_choose{1, 5}(i)
                if centernode_n ~= i
                    synchronize_pair{end + 1} = [centernode_n, i];
                    da_seg1_n(i) = NaN;  % 如果两振子同步了，就不再认为二者之间存在连边了
                else
                    da_seg1_n(i) = NaN;  % 如果两振子同步了，就不再认为二者之间存在连边了
                end
            end
        end
    end

    % % 第一步-稳定性筛选：计算第centernode_n个利用segment_num段信号计算的(dp, da)重叠图上同序号点到三点质心的最大欧氏距离max_distance_centroid
    % % 小于md则认为分支节点接收到的中心节点的平均模式是稳定的，因此认为二者之间存在连边
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
    
    % 新增：第一步-稳定性筛选：计算第centernode_n个利用segment_num段信号计算的(dp,da)重叠图上同序号点到三点质心的距离的方差
    % 如果该方差小于mcr则认为分支节点接收到的中心节点的平均模式是稳定的，因此认为二者之间存在连边
    centroid_rms = zeros(1, N);
    for i = 1:N
        points = [phase_difference_choose{1, centernode_n}(i), da_seg1_n(i)];
        for s = 2:segment_num
            points = [points; [phase_difference_choose{s, centernode_n}(i), amplitude_difference_choose{s, centernode_n}(i)]];
        end
       centroid_rms(i) = rmsOfDistanceFromCentroid(points); 
    end
    
    % 新增：保存centroid_rms从小到大排列后原本对应的序号
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

    %↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓如果只进行稳定性筛选，注释此段↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
    % 第二步-幅度筛选：基于幅度衰减信息进一步用聚类算法筛选出衰减较小的节点
    current_cluster = [];
    current_cluster_size = 0;
    for j = 1:N
        former_size = length(current_cluster);
        for i = 1:N
            current_data = da_seg1_n(i);
            if isempty(current_cluster)
                current_cluster = [max(da_seg1_n)];
            else
                % 计算当前数据与当前聚类中心的距离
                distance = abs(current_data - mean(current_cluster));
                % 如果距离小于等于聚类距离，将数据加入当前聚类
                if distance <= acd && distance ~= 0
                    current_cluster = [current_cluster, current_data];
                end
            end
        end
        %在数据量较少的情况下，需要补充扫描一次防止遗漏
        current_cluster = unique(current_cluster);
        current_cluster_size = length(current_cluster);
        if current_cluster_size == former_size
            break
        end
    end
    % 根据以上两步的结果构造邻接矩阵的第centernode_n行
    inference_adjacency_matrix(centernode_n, :) = ismember(da_seg1_n, current_cluster);
    inference_adjacency_matrix(centernode_n, centernode_n) = 0;
    %↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑如果只进行稳定性筛选，注释此段↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑
end
inference_adjacency_matrix(isnan(inference_adjacency_matrix)) = 0;

% 第三步-对称性筛选：由于MS方法检验的是某一节点的运动对剩余节点动力学行为的影响，因此结果是单向的
% 如果网络是无向网络，则将矩阵中所有的单向边双向化
for i = 1:N
    for j = 1:N
        if inference_adjacency_matrix(i, j) ~= inference_adjacency_matrix(j, i)  % 如果推断矩阵不对称
            %保留
            % inference_adjacency_matrix(i, j) = 1;
            % inference_adjacency_matrix(j, i) = 1;
            %剔除
            inference_adjacency_matrix(i, j) = 0;
            inference_adjacency_matrix(j, i) = 0;
        end
    end
end

% 第四步-为了适应耦合强度推断算法的需要，默认同步的节点之间存在连边
% 同时，在节点层面要将同步节点的所有连边共享
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

%计算推断的准确性
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
beta = 1;  %F-score的beta值，>1表示更看重recall，也就是希望漏判的目标边(FN)更少
F_beta_score = ((1+beta^2)*TP)/((1+beta^2)*TP+beta^2*FN+FP)
recall_node = TP_node ./ (TP_node + FN_node);
precision_node = TP_node ./ (TP_node + FP_node);
F_beta_score_node = ((1+beta^2)*TP_node)./((1+beta^2)*TP_node+beta^2*FN_node+FP_node);

% 根据邻接矩阵绘制推断目标网络
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
f1 = plot(G1, 'Layout', 'circle', 'NodeCData', nColors1); % 'force'布局可根据节点之间的连接力来布局节点
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
node_index = 1:N;
label = arrayfun(@num2str, node_index, 'UniformOutput', false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colormap parula
clim(common_color_limits)
labelnode(f1, 1:N, label)
title('推断目标网络')
f1.NodeFontWeight = 'bold';
f1.NodeLabelColor = 'r';
colorbar
subplot(1, 2, 2)
f2 = plot(G2, 'Layout', 'circle', 'NodeCData', nColors2); % 'force'布局可根据节点之间的连接力来布局节点
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_beta_score_label = arrayfun(@num2str, round(F_beta_score_node, 2), 'UniformOutput', false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colormap parula
clim(common_color_limits)
labelnode(f2, 1:N, F_beta_score_label)
title('推断结果网络')
f2.NodeFontWeight = 'bold';
f2.NodeLabelColor = 'r';
colorbar
% 双向图注释下面这部分
% f2.ArrowSize = 15;
% f2.EdgeLabel = G2.Edges.Weight;
% f2.EdgeFontSize = 10;
% f2.EdgeLabelColor = 'black';

%%
%%%%%%%%%%合并同步团后重新计算推断的准确性%%%%%%%%%%
% 新增：计算同步团内信号的PLV
synchronize_frequency_diff = [];
synchronize_PLV = [];
for i = 1:size(synchronize_pair, 2)
    synchronize_frequency_pair{i} = Average_frequency(synchronize_pair{i}, 1);
    synchronize_frequency_diff = [synchronize_frequency_diff, ...
        2*abs(synchronize_frequency_pair{i}(1)-synchronize_frequency_pair{i}(2))/ ...
        (synchronize_frequency_pair{i}(1)+synchronize_frequency_pair{i}(2))];
    synchronize_PLV = [synchronize_PLV, PLV(CL_Trajectory(:, synchronize_pair{i}(1)), ...
        CL_Trajectory(:, synchronize_pair{i}(2)))];  % 新增：计算同步对的PLV
end
Average_synchronize_frequency_diff = mean(synchronize_frequency_diff);
Average_synchronize_PLV = mean(synchronize_PLV);

% 新增：计算所有信号之间PLV
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

% 新增：根据频率去除错误判断的同步对
keep_index = [];
for sfp = synchronize_frequency_pair
    if abs(sfp{1}(1)-sfp{1}(2)) > 1 % 添加一层双保险
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

%计算推断的准确性
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
beta = 1;  %F-score的beta值，>1表示更看重recall，也就是希望漏判的目标边(FN)更少
F_beta_score_m = ((1+beta^2)*TP_m)/((1+beta^2)*TP_m+beta^2*FN_m+FP_m)
recall_m_node = TP_m_node ./ (TP_m_node + FN_m_node);
precision_m_node = TP_m_node ./ (TP_m_node + FP_m_node);
F_beta_score_m_node = ((1+beta^2)*TP_m_node)./((1+beta^2)*TP_m_node+beta^2*FN_m_node+FP_m_node);

% 根据同步邻接矩阵绘制推断目标网络
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
f3 = plot(G3, 'Layout', 'circle', 'NodeCData', nColors3); % 'force'布局可根据节点之间的连接力来布局节点
colormap parula
clim(common_color_limits)
labelnode(f3, 1:merge_N, mergelabel)
title('Target network after merging synchronize group', Interpreter='latex')
f3.NodeFontWeight = 'bold';
f3.NodeLabelColor = 'r';
colorbar
% 设置坐标轴标签字体样式
set(gca,'fontsize', 12);
subplot(1, 2, 2)
f4 = plot(G4, 'Layout', 'circle', 'NodeCData', nColors4); % 'force'布局可根据节点之间的连接力来布局节点
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_beta_score_m_label = arrayfun(@num2str, round(F_beta_score_m_node, 2), 'UniformOutput', false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colormap parula
clim(common_color_limits)
labelnode(f4, 1:merge_N, F_beta_score_m_label)
title('Result network after merging synchronize group', Interpreter='latex')
f4.NodeFontWeight = 'bold';
f4.NodeLabelColor = 'r';
colorbar
% 设置坐标轴标签字体样式
set(gca,'fontsize', 12);

toc

function Rms = rmsOfDistanceFromCentroid(points)
    % points 是一个包含n行2列的矩阵，每一行代表一个二维点的坐标 (x, y)
    % 计算质心
    centroid = mean(points, 1); % 按列求平均得到质心坐标
    % 计算每个点到质心的距离
    distances = pdist2(points, centroid);
    % 找出距离质心最远的点的距离
    Rms = rms(distances);
end

function synGroup = synPair2synGroup(synPair)
    % inputGroups: 输入的数组组成的单元格数组
    synGroup = {};

    while ~isempty(synPair)
        currentGroup = synPair{1};
        synPair(1) = [];

        % 查找与当前组有交集的其他组
        intersectingIndices = findIntersectingGroups(currentGroup, synPair);

        % 合并有交集的组
        while ~isempty(intersectingIndices)
            currentGroup = unique([currentGroup, cell2mat(synPair(intersectingIndices))]);
            synPair(intersectingIndices) = [];
            intersectingIndices = findIntersectingGroups(currentGroup, synPair);
        end

        % 将当前组添加到结果中
        synGroup{end + 1} = currentGroup;
    end
end

function mergeLabel = generateMergeLabel(node_index, nodeGroupsToMerge)
    % 生成同步团合并后的节点标签
    % 合并node_index和nodeGroupsToMerge并剔除掉nodeGroupsToMerge中包含的node_index的元素

    % 初始化结果数组
    mergeLabel = cell(1, length(nodeGroupsToMerge) + length(node_index));

    % 处理label的元素
    for i = 1:length(node_index)
        if ~any(cellfun(@(x) ismember(node_index(i), x), nodeGroupsToMerge))
            mergeLabel{i} = num2str(node_index(i));
        end
    end

    % 处理nodeGroupsToMerge的元素
    for i = 1:length(nodeGroupsToMerge)
        str = ['{', strjoin(arrayfun(@num2str, nodeGroupsToMerge{i}, 'UniformOutput', false), ', '), '}'];
        mergeLabel{length(node_index) + i} = str;
    end

    % 移除空的元素
    mergeLabel = mergeLabel(~cellfun('isempty', mergeLabel));
end

function mergedAdjMatrix = mergeSynchronizedNodeGroups(originalAdjMatrix, nodeGroupsToMerge, uniOption)
    % originalAdjMatrix: 原始的邻接矩阵
    % nodeGroupsToMerge: 要合并的节点组，每行是一个节点组
    % uniOption：逻辑值，true表示合并后的邻接矩阵只表示边的存在，元素只包含{0,1}
    % ..................false表示合并后的邻接矩阵中的耦合强度为原始矩阵对应边上耦合强度的和

    % 复制原始邻接矩阵
    mergedAdjMatrix = originalAdjMatrix;

    % 遍历每个节点组，合并节点并保留连边信息
    for i = 1:size(nodeGroupsToMerge, 2)
        nodeGroup = nodeGroupsToMerge{i};

        % 合并节点
        mergedNode_row = sum(mergedAdjMatrix(nodeGroup, :), 1);
        mergedNode_column = sum(mergedAdjMatrix(:, nodeGroup), 2);
        
        if uniOption == true
            % 将非零值设置为 1
            mergedNode_row = mergedNode_row ~= 0;
            mergedNode_column = mergedNode_column ~= 0;
        end

        % 将合并后的节点添加到邻接矩阵
        mergedAdjMatrix = [mergedAdjMatrix, mergedNode_column];
        mergedAdjMatrix = [mergedAdjMatrix; [mergedNode_row, 0]];

        % 移除原始节点
        mergedAdjMatrix(nodeGroup, :) = [];
        mergedAdjMatrix(:, nodeGroup) = [];

        % 更新未合并的同步团中的节点序号
        for j = i+1:size(nodeGroupsToMerge, 2)
            nodeGroupsToMerge{j} = updateNodeIndices(nodeGroupsToMerge{j}, nodeGroup);
        end
    end
end

function intersectingIndices = findIntersectingGroups(group, groups)
    % 查找与给定组有交集的其他组的索引
    intersectingIndices = [];
    for i = 1:length(groups)
        if any(ismember(groups{i}, group))
            intersectingIndices = [intersectingIndices, i];
        end
    end
end

function summedMatrix = sumMatrixByStep(matrix, jump_step)
    % 确定新矩阵的尺寸
    numRows = ceil(size(matrix, 1) / jump_step);
    newMatrix = zeros(numRows, size(matrix, 2));
    
    % 循环遍历原矩阵，按jump_step进行求和
    for i = 1:numRows
        startRow = (i - 1) * jump_step + 1;
        endRow = min(i * jump_step, size(matrix, 1));
        newMatrix(i,:) = sum(matrix(startRow:endRow,:), 1);
    end
    
    summedMatrix = newMatrix;
end

function updatedGroup = updateNodeIndices(group, removedNodes)
    % 更新未合并的同步团中的节点序号

    % 获取未合并的同步团中的节点在前一同步团合并后，旧序号到新序号的映射
    map = zeros(1, length(group));
    for i = 1:length(group)
        map(i) = sum(removedNodes < group(i));
    end
    
    % 更新未合并的同步团中的节点序号
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
    %params中包含{N, a, r, b, Laplace_matrix}
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
