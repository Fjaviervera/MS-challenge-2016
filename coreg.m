

spm('defaults','FMRI');
spm_jobman('initcfg');

matlabbatch{1}.spm.spatial.coreg.estimate.ref = {ref};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {source};
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'ncc';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

inputs = cell{0,0}
spm_jobman('run', matlabbatch,inputs);

