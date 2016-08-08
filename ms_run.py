import argparse
import lesions_seg_functions as lsf
import os
from os import listdir,makedirs
from os.path import isfile, join, isdir,exists,splitext
import time
import shutil
import platform
parser = argparse.ArgumentParser(description='Laimbio MS Challenge 2016 Method ')

parser.add_argument("t1_raw", help="path of the T1-w volume")
parser.add_argument("flair_raw", help="path of the FLAIR volume")
parser.add_argument("t1_pp", help="path of the T1-w volume")
parser.add_argument("flair_pp", help="path of the FLAIR volume")
parser.add_argument("brain_mask", help="path of the brain-mask volume")
parser.add_argument("spm_path",help="path to spm.sh if linux or spm.exe if windows")
parser.add_argument("output_folder", help="path where outputs and intermediate files will be stored")
parser.add_argument("-mcr_path", help="matlab compiler runtime path ")
parser.add_argument("-t2_raw", help="path of the T2-w volume")
parser.add_argument("-dp_raw", help="path of the DP volume")
parser.add_argument("-gado_raw", help="path of the T1-Gado volume")
parser.add_argument("-t2_pp", help="path of the T2-w volume")
parser.add_argument("-dp_pp", help="path of the DP volume")
parser.add_argument("-gado_pp", help="path of the T1-Gado volume")





args = parser.parse_args()
unpp_list = []
pp_list = []
if args.t2_raw == None or args.dp_raw == None or args.gado_raw == None or args.t2_pp == None or args.dp_pp == None or args.gado_pp == None:

        mode = 'basic'


        unpp_list.append(args.t1_raw)
        unpp_list.append(args.flair_raw)


        pp_list.append(args.t1_pp)
        pp_list.append(args.flair_pp)

        types_sequences= ['T1','FLAIR']
elif  args.dp_raw == None or args.gado_raw == None or args.dp_pp == None or args.gado_pp == None:

        mode = 'basic_t2'


        unpp_list.append(args.t1_raw)
        unpp_list.append(args.t2_raw)
        unpp_list.append(args.flair_raw)


        pp_list.append(args.t1_pp)
        unpp_list.append(args.t2_pp)
        pp_list.append(args.flair_pp)

        types_sequences= ['T1','T2','FLAIR']


elif  args.gado_raw == None or args.gado_pp == None:

    mode = 'basic_t2_dp'

    unpp_list.append(args.t1_raw)
    unpp_list.append(args.t2_raw)
    unpp_list.append(args.flair_raw)
    unpp_list.append(args.dp_raw)
    pp_list.append(args.t1_pp)
    unpp_list.append(args.t2_pp)
    pp_list.append(args.flair_pp)

    types_sequences = ['T1','T2','FLAIR','DP']

else:
        mode = 'all'


        unpp_list.append(args.t1_raw)
        unpp_list.append(args.t2_raw)
        unpp_list.append(args.flair_raw)


        unpp_list.append(args.dp_raw)
        unpp_list.append(args.gado_raw)

        pp_list.append(args.t1_pp)
        pp_list.append(args.t2_pp)
        pp_list.append(args.flair_pp)


        pp_list.append(args.dp_pp)
        pp_list.append(args.gado_pp)

        types_sequences = ['T1','T2', 'FLAIR' ,'DP','GADO']
mode_print =''
for seq in types_sequences:
    mode_print+=' '+ seq


print 'Starting MS lesion segmentation using '+ mode_print
print 'OS '+ platform.system()

subject=lsf.subject_class(types_sequences,unpp_list,args.output_folder)

if not exists(
        join(subject.dir, 'intermediate')):
    makedirs(
        join(subject.dir, 'intermediate'))

subject.add_intermediate(join(subject.dir, 'intermediate'))
subject.add_pp_sequences( types_sequences,pp_list)

if args.brain_mask == None:

    if args.t1_pp == None:
        subject = lsf.gunzip_T1_unpp(subject)

        if args.t1_raw[-3::] == '.gz':
            subject = lsf.gunzip_T1_unpp(subject)


        with open(join(subject.intermediate_path, "skull_batch.m"), "w") as f1:
            f1.write(' path_t1=\'' + subject.T1_gunzip_path + '\' ; \n ')
            f1.write('path_tpm=\'' + join(os.path.dirname(args.spm_path ),'spm12_mcr/spm12/tpm/TPM.nii') + '\' ; \n')
            with open(join(os.path.dirname(os.path.realpath(__file__)), "skullseg_creator.m")) as f:
                for line in f:
                    f1.write(line)

        if platform.system() == 'Windows':

            os.system(join(args.spm_path) + ' batch ' + join(subject.intermediate_path, "skull_batch.m"))

        else:
            if args.mcr_path == None:
                print ' Matlab runtime compiler path is needed to run in Linux/MacOS systems'
            os.system(args.spm_path + ' ' + args.mcr_path + ' batch ' + join(subject.intermediate_path, "skull_batch.m"))

    lsf.create_brain_mask(subject)

else:
    subject.add_brain_mask(args.brain_mask)

if args.t1_pp[-3::] == '.gz':
    subject=lsf.gunzip_T1_pp(subject)



shutil.copy(subject.T1_pp_gunzip_path,splitext(subject.T1_pp_gunzip_path)[0]+'_copia.nii')
subject.T1_pp_gunzip_path_copy = splitext(subject.T1_pp_gunzip_path)[0]+'_copia.nii'

with open(join(subject.intermediate_path, "coreg_batch.m"), "w") as f1:
    f1.write('ref=\'' + join(os.path.dirname(args.spm_path), 'spm12_mcr/spm12/canonical/avg305T1.nii') + '\' ; \n')
    f1.write('source=\'' + subject.T1_pp_gunzip_path + '\' ; \n ')
    with open(join(os.path.dirname(os.path.realpath(__file__)), "coreg.m")) as f:
        for line in f:
            f1.write(line)

if platform.system() == 'Windows':

    os.system(join(args.spm_path) + ' batch ' + join(subject.intermediate_path, "coreg_batch.m"))

else:
    if args.mcr_path == None:
        print ' Matlab runtime compiler path is needed to run in Linux/MacOS systems'
    os.system(args.spm_path + ' ' + args.mcr_path + ' batch ' + join(subject.intermediate_path, "coreg_batch.m"))


with open(join(subject.intermediate_path,"tissue_batch.m"), "w") as f1:
    f1.write(' path_t1=\'' + subject.T1_pp_gunzip_path+'\' ; \n ')
    f1.write('path_tpm=\'' + join(os.path.dirname(args.spm_path ),'spm12_mcr/spm12/tpm/TPM.nii')+'\' ; \n')
    with open(join(os.path.dirname(os.path.realpath(__file__)) ,"tissueseg_creator.m")) as f:

            for line in f:
                f1.write(line)

if platform.system() =='Windows':

    os.system(join(args.spm_path)+ ' batch '+ join(subject.intermediate_path,"tissue_batch.m" ))

else:
    if args.mcr_path ==None:
        print ' Matlab runtime compiler path is needed to run in Linux/MacOS systems'
    os.system(args.spm_path +' '+ args.mcr_path + ' batch ' + join(subject.intermediate_path, "tissue_batch.m"))


subject.add_tissue_segmentation( join(subject.intermediate_path,'c1T1_pp.nii'),join(subject.intermediate_path,'c2T1_pp.nii'),join(subject.intermediate_path,'c3T1_pp.nii'))


with open(join(subject.intermediate_path, "coreg_return_batch.m"), "w") as f1:
    f1.write(' ref=\'' + subject.T1_pp_gunzip_path_copy + '\' ; \n ')
    f1.write(' source=\'' + subject.T1_pp_gunzip_path + '\' ; \n ')
    f1.write(' tissue=\'' + subject.GM_path + '\' \n ' + '\''+subject.WM_path + '\' \n' + '\''+subject.CSF_path + '\' \n ;' )
    with open(join(os.path.dirname(os.path.realpath(__file__)), "coreg_return.m")) as f:
        for line in f:
            f1.write(line)

if platform.system() == 'Windows':

    os.system(join(args.spm_path) + ' batch ' + join(subject.intermediate_path, "coreg_return_batch.m"))

else:
    if args.mcr_path == None:
        print ' Matlab runtime compiler path is needed to run in Linux/MacOS systems'
    os.system(args.spm_path + ' ' + args.mcr_path + ' batch ' + join(subject.intermediate_path, "coreg_return_batch.m"))




subject=lsf.intensity_correction(subject)

subject=lsf.create_features_csfext(subject,2)


class_ext_csf_path= join(os.path.dirname(os.path.realpath(__file__)), 'classifiers','classifiers_server_csfext','classifier_200_None_'+mode,'clasifier_csf_ext_200_None.pkl')
subject.add_ext_csf_classifier(class_ext_csf_path)


subject=lsf.test_csfext(subject,n_estim=200, depth=None,flag=2)


subject=lsf.create_features_ms(subject,2)

class_ms_path= join(os.path.dirname(os.path.realpath(__file__)), 'classifiers','classifiers_server_ms','classifier_200_None_'+mode,'clasifier_ms_200_None.pkl')
subject.add_ms_classifier(class_ms_path)


subject=lsf.test_ms(subject,n_estim=200, depth=None,flag=2)
if mode == 'basic':
    subject=lsf.lesion_growing(subject,theta=0.15,beta_grow=2,flag=1)
else:
    subject = lsf.lesion_growing(subject, theta=0.25, beta_grow=2, flag=1)

print 'Job done  '
