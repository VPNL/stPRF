classdef Constants
    % 	properties( Constant = true )
    %
    %     end
    %
    
    methods (Static)
        function d = getDir
            
            output_dir = '/share/kalanit/users/insubkim/Documents/demo_test/stPRF/results/example/data';
            simulator_dir = '/share/kalanit/users/insubkim/Documents/demo_test/stPRF/results/example/data';
            sessionDate = 'simulation';
            d.analysisoption = 9; 
            
            d.sessionDate = sessionDate;
            d.output_dir = output_dir;
            d.simulator_dir = simulator_dir;
            d.IRF_dir = fullfile(output_dir,'IRF',sessionDate,'/');
            d.grid_dir = fullfile(output_dir,'grid',sessionDate,'/');
            d.plot_dir = fullfile(output_dir,'plot',sessionDate,'/');
            

            % simulaiton related
            [~,id] = fileparts(fileparts(simulator_dir));
            d.sim_json_dir = fullfile(simulator_dir,'jsonfiles',sessionDate,'/');
            d.sim_grid_dir = fullfile(simulator_dir,'grid',sessionDate,'/');
            d.sim_IRF_dir = fullfile(simulator_dir,'IRF',sessionDate,'/');
            d.sim_BOLD_dir = fullfile(simulator_dir,'synBOLD',sessionDate,'/');
            d.sim_pt = fullfile(simulator_dir,'dataTable',sessionDate,'/','gray',['df_' id '.mat']);

            
        end
        
        function s = getSub
            nativeOakProjectPth = '/oak/stanford/groups/kalanit/biac2/kgs/projects';
            anatNames  = {'kgs201310_v53','kis202008_v60','ek202011_v60', ...
                'hk2019_v60', 'lm202105_v60', 'ac202105_v60'...
                'jj202105_v60','jr202105_v60','xy202105_v60', 'yc202106_v60'};
            usedSessions = {'session1', 'session2', 'session1', ...
                'session1', 'session1', 'session1', ...
                'session1', 'session1', 'session1', 'session1'};
            
            for SN = 1:length(anatNames)
                subjID   = sprintf('subj%02d', SN);
                anatName = anatNames{SN}; session = usedSessions{SN};
                
                sessionDir{SN} = fullfile(nativeOakProjectPth, ...
                    './spatiotemporal','experiments','stRet',[ 'data/' subjID '/' session]);
            end
            s.sessionDir = sessionDir;

        end
        
    end
    
    
    
end