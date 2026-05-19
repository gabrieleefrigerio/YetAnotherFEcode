classdef FrequencyAnalyzer < handle
    % FREQUENCYANALYZER Classe per l'analisi in frequenza (FFT e FRF) 
    % delle risposte dinamiche transitorie.
    
    properties
        sim_data        % Riferimento all'oggetto TransientSimulator
        A_mag           % Ampiezza dell'impulso [N]
        duration_pulse  % Durata dell'impulso [s]
        
        % Assi
        Fs              % Frequenza di campionamento [Hz]
        L               % Numero di campioni
        f               % Vettore delle frequenze [Hz]
        time_vec        % Vettore temporale [s]
        
        % Segnali nel tempo
        F_input         % Storia temporale della forza ricostruita
        
        % Spettri (Ampiezze Single-Sided)
        P1_F            % Spettro della Forza
        P1_X_fom        % Spettro dello Spostamento (FOM)
        P1_X_rom        % Spettro dello Spostamento (ROM)
        
        % Funzioni di Risposta in Frequenza (Magnitudo in dB)
        FRF_fom_dB
        FRF_rom_dB
    end
    
    methods
        function obj = FrequencyAnalyzer(transient_simulator, A_mag, duration_pulse)
            % Costruttore: Inizializza l'oggetto ed esegue subito i calcoli
            obj.sim_data = transient_simulator;
            obj.A_mag = A_mag;
            obj.duration_pulse = duration_pulse;
            
            % Avvia l'elaborazione dei dati
            obj.process_data();
        end
        
        function process_data(obj)
            % 1. Setup assi tempi e frequenze
            obj.time_vec = obj.sim_data.time_fom;
            dt = obj.time_vec(2) - obj.time_vec(1); % Passo temporale (h)
            obj.Fs = 1 / dt;
            obj.L = length(obj.time_vec);
            obj.f = obj.Fs * (0:(obj.L/2)) / obj.L;
            
            % 2. Ricostruzione Forza (Mezzo seno)
            obj.F_input = zeros(size(obj.time_vec));
            idx_pulse = obj.time_vec <= obj.duration_pulse;
            obj.F_input(idx_pulse) = obj.A_mag * sin(pi * obj.time_vec(idx_pulse) / obj.duration_pulse);
            
            % 3. Calcolo Spettri (FFT) tramite metodo interno
            obj.P1_F     = obj.compute_single_sided_fft(obj.F_input);
            obj.P1_X_fom = obj.compute_single_sided_fft(obj.sim_data.disp_center_fom);
            obj.P1_X_rom = obj.compute_single_sided_fft(obj.sim_data.disp_center_rom);
            
            % 4. Calcolo FRF Empirica in Decibel (dB)
            % Aggiungiamo 'eps' per evitare divisioni per zero
            FRF_fom = obj.P1_X_fom ./ (obj.P1_F + eps);
            FRF_rom = obj.P1_X_rom ./ (obj.P1_F + eps);
            
            obj.FRF_fom_dB = 20 * log10(FRF_fom + eps);
            obj.FRF_rom_dB = 20 * log10(FRF_rom + eps);
        end
        
        function P1 = compute_single_sided_fft(obj, signal)
            % UTILITY: Calcola lo spettro normalizzato (Single-Sided) di un segnale
            Y = fft(signal);
            P2 = abs(Y / obj.L);
            P1 = P2(1:floor(obj.L/2)+1);
            P1(2:end-1) = 2 * P1(2:end-1);
            
            % Forza l'output come vettore colonna per coerenza nei plot
            if isrow(P1)
                P1 = P1'; 
            end
        end
        
        function plot_results(obj, f_max_plot)
            % PLOT_RESULTS Genera la figura a 3 pannelli con spettri e FRF.
            % f_max_plot: (Opzionale) Frequenza massima da visualizzare sull'asse X
            
            if nargin < 2
                f_max_plot = 1500; % Default: visualizza fino a 1500 Hz
            end
            
            % Forza f come vettore colonna per il plot
            f_plot = obj.f;
            if isrow(f_plot), f_plot = f_plot'; end
            
            figure('Name', 'Frequency Domain Analysis', 'Color', 'w', 'Position', [150, 100, 800, 800]);
            
            % --- Subplot 1: Input ---
            subplot(3, 1, 1);
            semilogy(f_plot, obj.P1_F, 'k-', 'LineWidth', 1.5);
            xlim([0, f_max_plot]);
            grid on;
            title('Impulse Spectrum (Input Force)', 'FontSize', 11);
            ylabel('Force [N]');
            xlabel('Frequency [Hz]');
            
            % --- Subplot 2: Output ---
            subplot(3, 1, 2);
            semilogy(f_plot, obj.P1_X_fom * 1000, 'b-', 'LineWidth', 1.5); hold on;
            semilogy(f_plot, obj.P1_X_rom * 1000, 'r--', 'LineWidth', 1.5);
            xlim([0, f_max_plot]);
            grid on;
            title('Displacement Spectrum (Output Responses)', 'FontSize', 11);
            ylabel('Displacement [mm]');
            xlabel('Frequency [Hz]');
            legend({'FOM', 'ROM[1 3 5 7 9 11 MC]'}, 'Location', 'northeast');
            
            % --- Subplot 3: FRF ---
            subplot(3, 1, 3);
            plot(f_plot, obj.FRF_fom_dB, 'b-', 'LineWidth', 1.5); hold on;
            plot(f_plot, obj.FRF_rom_dB, 'r--', 'LineWidth', 1.5);
            xlim([0, f_max_plot]);
            % Taglia l'asse Y in basso per escludere il rumore numerico
            ylim([min(obj.FRF_fom_dB), max(obj.FRF_fom_dB)+10]); 
            grid on;
            title('Frequency Response Function (FRF Magnitude)', 'FontSize', 11);
            ylabel('Magnitude [dB]');
            xlabel('Frequency [Hz]');
            legend({'FOM', 'ROM[1 3 5 7 9 11 MC]'}, 'Location', 'northeast');
            
            sgtitle('Spectral Analysis: FOM vs ROM[1 3 5 7 9 11 MC]', 'FontSize', 14, 'FontWeight', 'bold');
        end
    end
end