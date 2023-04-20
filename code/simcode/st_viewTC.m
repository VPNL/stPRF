function st_viewTC(DT,voxelNumber)
figure('position',[500 500 1200 600])

% CST model
subplot(321);
plot(DT{1}(voxelNumber,:).tc)
model =DT{1}(voxelNumber,:).tmodel{1};
title(sprintf('%s noiseless',model));
xlabel('time'); ylabel('amp');

subplot(322);
plot(DT{1}(voxelNumber,:).ntc)
title(sprintf('%s + noise',model))
xlabel('time'); ylabel('amp');

% DN model
subplot(323);
plot(DT{2}(voxelNumber,:).tc)
model=DT{2}(voxelNumber,:).tmodel{1};
title(sprintf('%s noiseless',model));
xlabel('time'); ylabel('amp');

subplot(324);
plot(DT{2}(voxelNumber,:).ntc)
title(sprintf('%s + noise',model))
xlabel('time'); ylabel('amp');

% spatial model
subplot(325);
plot(DT{3}(voxelNumber,:).tc)
model=DT{3}(voxelNumber,:).tmodel{1};
title(sprintf('%s noiseless',model));
xlabel('time'); ylabel('amp');

subplot(326);
plot(DT{3}(voxelNumber,:).ntc)
title(sprintf('%s + noise',model))
xlabel('time'); ylabel('amp');

disp([DT{1}(voxelNumber,:); DT{2}(voxelNumber,:); DT{3}(voxelNumber,:)])

end