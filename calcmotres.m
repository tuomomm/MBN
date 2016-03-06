%calcmotres.m: Calculate the single motif promotion figure.
%Tuomo M?ki-Marttunen, 2013-2016

Nsamp = 20;

N = 100;
ps = 0.025:0.025:0.5;
Mots = zeros(Nsamp,16,20,16);
RandMots = zeros(Nsamp,16,20,16);
load(['Mots_mean.mat'],'Motsm','RandMotsm','Motsd','RandMotsd','Motsmp','Motsdp','RandMotsmp','RandMotsdp');
tic;

for ITER = 1:Nsamp
  rand('state',ITER);randn('state',ITER);
  for ip=1:length(ps)
    for imot=1:16
      M = mbn(N,ps(ip),[zeros(1,imot-1) 1 zeros(1,16-imot)],1,inf,1);
      Mots(ITER,:,ip,imot) = givemotifs(M);
      M = rand(100) <= ps(ip); M = M-diag(diag(M));
      RandMots(ITER,:,ip,imot) = givemotifs(M);
    end
    disp(['ITER = ' num2str(ITER) ', p = ' num2str(ps(ip)) ' done']);
  end
  toc
end

Motsm = nanmean(Mots);
Motsd = nanstd(Mots);
Motsmp = nanmean(Mots,3);
Motsdp = nanstd(Mots,[],3);
RandMotsm = nanmean(RandMots);
RandMotsd = nanstd(RandMots);
RandMotsmp = nanmean(RandMots,3);
RandMotsdp = nanstd(RandMots,[],3);
save(['Mots_mean.mat'],'Motsm','RandMotsm','Motsd','RandMotsd','Motsmp','Motsdp','RandMotsmp','RandMotsdp');


