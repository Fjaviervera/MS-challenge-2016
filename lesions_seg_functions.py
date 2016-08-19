from medpy.io import load
from medpy.io import save
import numpy as np
from sklearn import utils
from os import listdir,makedirs
from os.path import isfile, join, isdir,exists
import os
from medpy.features import indices
import pickle
import medpy.metric
import scipy.ndimage as ndimage
import math
from sklearn.externals import joblib
from sklearn.ensemble import RandomForestClassifier
from sklearn import mixture
import scipy.stats as stats
from medpy import filter
import matplotlib.pyplot as plt
import xlsxwriter
from scipy.ndimage.morphology import  binary_dilation





class subject_class:


    def __init__(self, types_sequences ,path_sequences,path_subject):


        self.dir = path_subject
        self.T1= 0
        self.T1_path = ''
        self.T2 = 0
        self.T2_path= ''
        self.FLAIR = 0
        self.FLAIR_path = ''
        self.DP = 0
        self.DP_path = ''
        self.GADO = 0
        self.GADO_path = ''

        for n_paths in range(len(path_sequences)):

            if types_sequences [n_paths] =='T1':
                self.T1 = 1
                self.T1_path = path_sequences[n_paths]

            if types_sequences [n_paths] =='T2':
                self.T2 = 1
                self.T2_path = path_sequences[0]

            if types_sequences [n_paths] =='FLAIR':
                self.FLAIR = 1
                self.FLAIR_path = path_sequences[n_paths]

            if types_sequences [n_paths] =='GADO':
                self.GADO = 1
                self.GADO_path = path_sequences[n_paths]

            if types_sequences [n_paths] =='DP':
                self.DP = 1
                self.DP_path = path_sequences[n_paths]



        self.T1_irs = 0
        self.T1_irs_path = ''
        self.T2_irs = 0
        self.T2_irs_path = ''
        self.FLAIR_irs = 0
        self.FLAIR_irs_path = ''
        self.DP_irs = 0
        self.DP_irs_path = ''
        self.GADO_irs = 0
        self.GADO_irs_path = ''

        self.brain_mask = 0
        self.brain_mask_path = ''

        self.lesion_mask = 0
        self.lesion_mask_path = ''
        self.csf_ext_train = 0
        self.csf_ext_train_path = ''

    def add_intermediate(self,path_intermediate):
        self.intermediate = 1
        self.intermediate_path = path_intermediate

    def add_brain_mask(self,path_brain_mask):
        self.brain_mask = 1
        self.brain_mask_path = path_brain_mask

    def add_pp_sequences(self,types_sequences ,path_pp_sequences):

        for n_paths in range(len(path_pp_sequences)):

            if types_sequences[n_paths] == 'T1':
                self.T1_pp = 1
                self.T1_pp_path = path_pp_sequences[n_paths]

            if types_sequences[n_paths] == 'T2':
                self.T2_pp = 1
                self.T2_pp_path = path_pp_sequences[0]

            if types_sequences[n_paths] == 'FLAIR':
                self.FLAIR_pp = 1
                self.FLAIR_pp_path = path_pp_sequences[n_paths]

            if types_sequences[n_paths] == 'GADO':
                self.GADO_pp = 1
                self.GADO_pp_path = path_pp_sequences[n_paths]

            if types_sequences[n_paths] == 'DP':
                self.DP_pp = 1
                self.DP_pp_path = path_pp_sequences[n_paths]

    def add_tissue_segmentation(self, path_GM,path_WM,path_CSF):
        self.GM = 1
        self.GM_path = path_GM
        self.WM = 1
        self.WM_path = path_WM
        self.CSF = 1
        self.CSF_path = path_CSF



    def add_lesion_mask(self,path_lesion_mask):
        self.brain_mask = 1
        self.brain_mask = path_lesion_mask

    def return_sequences(self):
        sequences = []
        if self.T1 == 1:
            sequences.append(self.T1_path)
        if self.T2 == 1:
            sequences.append(self.T2_path)
        if self.FLAIR == 1:
            sequences.append(self.FLAIR_path)
        if self.DP == 1:
            sequences.append(self.DP_path)
        if self.GADO == 1:
            sequences.append(self.GADO_path)
        return sequences

    def return_irs_sequences(self):
        sequences = []
        if self.T1_irs == 1:
            sequences.append(self.T1_irs_path)
        if self.T2_irs == 1:
            sequences.append(self.T2_irs_path)
        if self.FLAIR_irs == 1:
            sequences.append(self.FLAIR_irs_path)
        if self.DP_irs == 1:
            sequences.append(self.DP_irs_path)
        if self.GADO_irs == 1:
            sequences.append(self.GADO_irs_path)
        return sequences

    def add_ext_csf_classifier(self,class_ext_csf_path):
        self.classifier_ext_csf = 1
        self.classifier_ext_csf_path = class_ext_csf_path

    def add_ms_classifier(self, class_ms_path):
        self.classifier_ms = 1
        self.classifier_ms_path = class_ms_path

def create_subject(flags_sequences,path_sequences,subject_name,subject_path):

    subject = subject_class(subject_name, flags_sequences, path_sequences,subject_path)

    return subject

def create_brain_mask(subject):
    if subject.T1_pp == 1:
        vol, header = load(subject.T1_pp_path)
        mask_brain = vol > 0
        save(mask_brain, join(subject.intermediate_path, 'brain_mask.nii.gz'), header)

    else:
        c1, headerc1 = load(join(subject.intermediate_path,'c1T1.nii'))
        c2, headerc2 = load(join(subject.intermediate_path, 'c2T1.nii'))
        c3, headerc3 = load(join(subject.intermediate_path, 'c3T1.nii'))
        c4, headerc4 = load(join(subject.intermediate_path, 'c4T1.nii'))
        c5, headerc5 = load(join(subject.intermediate_path, 'c5T1.nii'))
        t1,t1_header = load(subject.T1_path)


        mask_brain = ((c1+c2+c3) > 0.95) * (c4 + c5 < 0.1 )
        save(mask_brain, join(subject.intermediate_path,'brain_mask.nii.gz'), t1_header)
    subject.add_brain_mask(join(subject.intermediate_path, 'brain_mask.nii.gz'))
    return subject
def retrain_intensity(image,seq_type):

    path = 'train_int'
    subdirectories = os.listdir(path)
    folders = []
    for dir in subdirectories:
        if isdir(join(path, dir)):
            folders.append(dir)
    images=[]
    if seq_type == 1:
        seq = 'T1_preprocessed.nii.gz'
    if seq_type == 2:
        seq = 'T2_preprocessed.nii.gz'
    if seq_type == 3:
        seq = 'FLAIR_preprocessed.nii.gz'
    if seq_type == 4:
        seq = 'DP_preprocessed.nii.gz'
    if seq_type == 5:
        seq = 'GADO_preprocessed.nii.gz'
    for subject in folders:
        im, im_header = medpy.io.load(join(path, subject, seq))
        mask, mask_header = medpy.io.load(join(path, subject, 'Mask_registered.nii.gz'))

        images.append(im[mask > 0])

    images.append(image)
    irs = medpy.filter.IntensityRangeStandardization()

    trained_model, transformed_images = irs.train_transform(images)

    return transformed_images[15]

def intensity_correction(subject):

        # if not exists(join(subject.dir, subject.irs_dir)):
        #     makedirs(join(subject.dir, subject.irs_dir))
        # if not exists(join(subject.dir, subject.irs_dir,subject.name)):
        #     makedirs(join(subject.dir, subject.irs_dir,subject.name))

        if subject.T1 == 1:

            with open(join(os.path.dirname(os.path.realpath(__file__)),'models/model_T1.pkl'), 'r') as f:
                irs = pickle.load(f)

            vol, header = load(subject.T1_pp_path)
            mask, mask_header = load(subject.brain_mask_path)
            spacing_indices = (1, 1, 1)

            vol_indices = indices(vol, spacing_indices, mask)

            try:
                intensities_corrected = irs.transform(vol[mask > 0])
            except:
                intensities_corrected= retrain_intensity(vol[mask > 0], 1)


            vol_irs = np.zeros(vol.shape)

            for n_voxel in range(len(vol[mask > 0])):
                vol_irs[vol_indices[n_voxel][0], vol_indices[n_voxel][1], vol_indices[n_voxel][2]] = \
                    intensities_corrected[n_voxel]

            save(vol_irs, join(subject.intermediate_path, 'T1_corrected.nii.gz'), header)

            subject.T1_irs = 1
            subject.T1_irs_path = join(subject.intermediate_path, 'T1_corrected.nii.gz')

        if subject.T2 == 1:

            with open(join(os.path.dirname(os.path.realpath(__file__)),'models/model_T2.pkl'), 'r') as f:
                irs = pickle.load(f)

            vol, header = load(subject.T2_pp_path)
            mask, mask_header = load(subject.brain_mask_path)
            spacing_indices = (1, 1, 1)

            vol_indices = indices(vol, spacing_indices, mask)

            try:
                intensities_corrected = irs.transform(vol[mask > 0])
            except:
                intensities_corrected = retrain_intensity(vol[mask > 0], 2)

            vol_irs = np.zeros(vol.shape)

            for n_voxel in range(len(vol[mask > 0])):
                vol_irs[vol_indices[n_voxel][0], vol_indices[n_voxel][1], vol_indices[n_voxel][2]] = \
                    intensities_corrected[n_voxel]

            save(vol_irs, join(subject.intermediate_path, 'T2_corrected.nii.gz'), header)

            subject.T2_irs = 1
            subject.T2_irs_path = join(subject.intermediate_path,'T2_corrected.nii.gz')

        if subject.FLAIR == 1:

            with open(join(os.path.dirname(os.path.realpath(__file__)),'models/model_FLAIR.pkl'), 'r') as f:
                irs = pickle.load(f)

            vol, header = load(subject.FLAIR_pp_path)
            mask, mask_header = load(subject.brain_mask_path)
            spacing_indices = (1, 1, 1)

            vol_indices = indices(vol, spacing_indices, mask)

            try:
                intensities_corrected = irs.transform(vol[mask > 0])
            except:
                intensities_corrected = retrain_intensity(vol[mask > 0], 3)

            vol_irs = np.zeros(vol.shape)

            for n_voxel in range(len(vol[mask > 0])):
                vol_irs[vol_indices[n_voxel][0], vol_indices[n_voxel][1], vol_indices[n_voxel][2]] = \
                    intensities_corrected[n_voxel]

            save(vol_irs, join(subject.intermediate_path,'FLAIR_corrected.nii.gz'), header)

            subject.FLAIR_irs = 1
            subject.FLAIR_irs_path = join(subject.intermediate_path, 'FLAIR_corrected.nii.gz')

        if subject.DP == 1:

            with open(join(os.path.dirname(os.path.realpath(__file__)),'models/model_DP.pkl'), 'r') as f:
                irs = pickle.load(f)

            vol, header = load(subject.DP_pp_path)
            mask, mask_header = load(subject.brain_mask_path)
            spacing_indices = (1, 1, 1)

            vol_indices = indices(vol, spacing_indices, mask)

            try:
                intensities_corrected = irs.transform(vol[mask > 0])
            except:
                intensities_corrected = retrain_intensity(vol[mask > 0], 4)

            vol_irs = np.zeros(vol.shape)

            for n_voxel in range(len(vol[mask > 0])):
                vol_irs[vol_indices[n_voxel][0], vol_indices[n_voxel][1], vol_indices[n_voxel][2]] = \
                    intensities_corrected[n_voxel]

            save(vol_irs, join(subject.intermediate_path, 'DP_corrected.nii.gz'), header)

            subject.DP_irs = 1
            subject.DP_irs_path = join(subject.intermediate_path, 'DP_corrected.nii.gz')

        if subject.GADO== 1:

            with open(join(os.path.dirname(os.path.realpath(__file__)),'models/model_GADO.pkl'), 'r') as f:
                irs = pickle.load(f)

            vol, header = load(subject.GADO_pp_path)
            mask, mask_header = load(subject.brain_mask_path)
            spacing_indices = (1, 1, 1)

            vol_indices = indices(vol, spacing_indices, mask)

            try:
                intensities_corrected = irs.transform(vol[mask > 0])
            except:
                intensities_corrected = retrain_intensity(vol[mask > 0], 5)

            vol_irs = np.zeros(vol.shape)

            for n_voxel in range(len(vol[mask > 0])):
                vol_irs[vol_indices[n_voxel][0], vol_indices[n_voxel][1], vol_indices[n_voxel][2]] = \
                    intensities_corrected[n_voxel]

            save(vol_irs, join(subject.intermediate_path, 'GADO_corrected.nii.gz'), header)

            subject.GADO_irs = 1
            subject.GADO_irs_path = join(subject.intermediate_path, 'GADO_corrected.nii.gz')

        return subject

def gray_matter_threshold(flair,GM_mask):


    GM=medpy.features.intensities(flair,GM_mask)





    n, bins, patch = plt.hist(GM, bins=200, label="Flair WM",histtype='stepfilled')

    maximo = np.max(n)
    paso=bins[1]-bins[0]
    pico=bins[np.argmax(n)]+paso/2



    test=n>(maximo/2)
    iter=0
    modo=0


    for j in test:
        if modo==0:
            if j == True:
                indice1=iter
                modo=1
        else:
            if j==False:
                indice2=iter-1
                break
        iter+=1


    fwhm=(bins[indice2]-bins[indice1])/2




    # plt.axvline(pico, color='r', linestyle='dashed', linewidth=2)


    # gamma_int=1
    gamma_int =0.2

    # plt.axvline(pico+gamma_int*fwhm, color='r', linestyle='dashed', linewidth=5,label = "Threshold")





    Tint=pico + gamma_int*fwhm

    return flair>Tint

def kernel_subs_creator(tam,spacing):

    kernel_shape = np.asarray([(int(round(tam/spacing[0]))),(int(round(tam/spacing[1]))),(int(round(tam/spacing[2])))])
    for shape_kernel in range(len(kernel_shape)):
        if kernel_shape[shape_kernel] % 2 ==0:
                kernel_shape[shape_kernel] +=1
    kernel = np.ones(kernel_shape.tolist())*-1/((kernel_shape[0] * kernel_shape[1] *kernel_shape[2]) -1)
    kernel_center = [ math.floor(elem/2.0) for elem in kernel.shape ]

    kernel[int(kernel_center[0]),int(kernel_center[1]),int(kernel_center[2])] = 1
    return kernel

def subsampling(path,sampling):

    random_data=utils.shuffle(np.load(path,'r'))
    subset=random_data[1:int(sampling*random_data.shape[0])][:].copy()

    return subset

def create_features_csfext(subject,flag):

    print "------------Creating csf ext features-----------------"

    mask, mask_header = load(subject.brain_mask_path)
    pv1, pv1_header = load(subject.GM_path)
    pv2, pv2_header = load(subject.WM_path)
    pv0, pv0_header = load(subject.CSF_path)
    flair, fl_header = load(subject.FLAIR_irs_path)

    if flag == 0:
        mask_class, mask_class_header = load(subject.gt_extcsf_path)
        mask_class = ((pv0 > pv1) * (pv0 > pv2) * mask_class == 0)
        data_name = 'vent'

    if flag == 1:

        mask_class,mask_class_header = load(subject.gt_extcsf_path)
        data_name = 'csf_ext'


    if flag == 2:

        mask_class = ((pv0 > pv1) * (pv0 > pv2))
        data_name = 'csf'

    if flag > 3:

        'error invalid flag'
        return -1



    mask_voxels = mask * mask_class

    iter = 0

    for vol_path in subject.return_irs_sequences():


        print 'Generating csf features ---< ' + '  ' + vol_path

        image_data, image_header = load(vol_path)

        spacing = medpy.io.header.get_pixel_spacing(image_header)

        intensities = medpy.features.intensity.intensities(image_data, mask_voxels)

        gaussian1 = medpy.features.intensity.local_mean_gauss(image_data, 3, spacing, mask_voxels)
        gaussian2 = medpy.features.intensity.local_mean_gauss(image_data, 5, spacing, mask_voxels)
        gaussian3 = medpy.features.intensity.local_mean_gauss(image_data, 7, spacing, mask_voxels)

        kernel1 = kernel_subs_creator(3, spacing)
        subs1 = medpy.features.intensity.intensities(ndimage.filters.convolve(image_data, kernel1), mask_voxels)
        kernel2 = kernel_subs_creator(5, spacing)
        subs2 = medpy.features.intensity.intensities(ndimage.filters.convolve(image_data, kernel2), mask_voxels)
        kernel3 = kernel_subs_creator(7, spacing)
        subs3 = medpy.features.intensity.intensities(ndimage.filters.convolve(image_data, kernel3), mask_voxels)

        if iter == 0:
            joined = medpy.features.utilities.join(intensities, gaussian1, gaussian2, gaussian3, subs1, subs2, subs3)
        else:
            joined = medpy.features.utilities.join(joined, intensities, gaussian1, gaussian2, gaussian3, subs1, subs2,
                                                   subs3)
        iter += 1
    print 'Generating tissue and distance features'
    spacing_indices = (1, 1, 1)

    spacing = medpy.io.header.get_pixel_spacing(fl_header)

    indices = medpy.features.indices(flair, spacing_indices, mask_voxels)

    distances_0 = medpy.features.intensity.centerdistance_xdminus1(flair, 0, spacing, mask_voxels)
    distances_1 = medpy.features.intensity.centerdistance_xdminus1(flair, 1, spacing, mask_voxels)
    distances_2 = medpy.features.intensity.centerdistance_xdminus1(flair, 2, spacing, mask_voxels)



    flair_unpp, fl_unpp_header = load(subject.FLAIR_path)

    skull = (flair_unpp * (mask==0))>10


    spacing = medpy.io.header.get_pixel_spacing(fl_header)



    dist_transform_vol = ndimage.distance_transform_edt(skull == 0, sampling=spacing)
    dist_transform_feature = medpy.features.intensity.intensities(dist_transform_vol, mask_voxels)


    intensities_pv0 = medpy.features.intensity.intensities(pv0, mask_voxels)
    intensities_pv1 = medpy.features.intensity.intensities(pv1, mask_voxels)
    intensities_pv2 = medpy.features.intensity.intensities(pv2, mask_voxels)

    gaussian1_pv0 = medpy.features.intensity.local_mean_gauss(pv0, 3, spacing, mask_voxels)
    gaussian2_pv0 = medpy.features.intensity.local_mean_gauss(pv0, 7, spacing, mask_voxels)
    gaussian3_pv0 = medpy.features.intensity.local_mean_gauss(pv0, 15, spacing, mask_voxels)

    gaussian1_pv1 = medpy.features.intensity.local_mean_gauss(pv1, 3, spacing, mask_voxels)
    gaussian2_pv1 = medpy.features.intensity.local_mean_gauss(pv1, 7, spacing, mask_voxels)
    gaussian3_pv1 = medpy.features.intensity.local_mean_gauss(pv1, 15, spacing, mask_voxels)

    gaussian1_pv2 = medpy.features.intensity.local_mean_gauss(pv2, 3, spacing, mask_voxels)
    gaussian2_pv2 = medpy.features.intensity.local_mean_gauss(pv2, 7, spacing, mask_voxels)
    gaussian3_pv2 = medpy.features.intensity.local_mean_gauss(pv2, 15, spacing, mask_voxels)

    joined = medpy.features.utilities.join(joined, distances_0, distances_1, distances_2, intensities_pv0, intensities_pv1,
                                           intensities_pv2, gaussian1_pv0, gaussian2_pv0, gaussian3_pv0, gaussian1_pv1
                                           , gaussian2_pv1, gaussian3_pv1, gaussian1_pv2, gaussian2_pv2, gaussian3_pv2,
                                            dist_transform_feature)




    if flag ==0:

        if not exists(join(subject.dir, 'descriptors')):
            makedirs(join(subject.dir, 'descriptors'))

        np.save(join(subject.dir, 'descriptors', data_name + '_descriptor'), joined)

        np.save(join(subject.dir, 'descriptors', data_name + '_indices'), indices)

        subject.class0_extcsf_descriptor_path = join(subject.dir,'descriptors',data_name+'_descriptor.npy')
        subject.class0_extcsf_indices_path = join(subject.dir,'descriptors',data_name+'_indices.npy')
    if flag == 1:

        if not exists(join(subject.dir, 'descriptors')):
            makedirs(join(subject.dir, 'descriptors'))

        np.save(join(subject.dir, 'descriptors', data_name + '_descriptor'), joined)

        np.save(join(subject.dir, 'descriptors', data_name + '_indices'), indices)

        subject.class1_extcsf_descriptor_path = join(subject.dir, 'descriptors',
                                                     data_name + '_descriptor.npy')
        subject.class1_extcsf_indices_path = join(subject.dir, 'descriptors', data_name + '_indices.npy')


    if flag == 2:

        if not exists(join(subject.intermediate_path, 'descriptors')):
            makedirs(join(subject.intermediate_path, 'descriptors'))

        np.save(join(subject.intermediate_path, 'descriptors', data_name + '_descriptor'), joined)

        np.save(join(subject.intermediate_path, 'descriptors', data_name + '_indices'), indices)

        subject.test_extcsf_descriptor_path = join(subject.intermediate_path, 'descriptors',
                                                     data_name + '_descriptor.npy')
        subject.test_extcsf_indices_path = join(subject.intermediate_path, 'descriptors', data_name + '_indices.npy')
    print 'csf features saved'
    return subject

def train_subject_csfext(subject,subject_list,n_estim,depth,jobs):

    lista_0=[]
    lista_1=[]
    for subject_to_train in subject_list:

        if subject_to_train.name!=subject.name:
            print subject_to_train.name


            lista_0.append(subsampling(subject_to_train.class0_extcsf_descriptor_path,0.1))
            lista_1.append(subsampling(subject_to_train.class1_extcsf_descriptor_path, 0.3))



    print 'concatenando'
    clase0=np.concatenate(lista_0,0)
    lista_0=None


    clase1=np.concatenate(lista_1,0)

    lista_1=None

    print clase0.shape
    print clase1.shape

    label0=np.ones([clase0.shape[0],],dtype=int)

    label1=2*np.ones((clase1.shape[0],),dtype=int)




    print 'subset creado comenzando entrenamiento'


    clf = RandomForestClassifier(n_estimators=n_estim, max_depth=depth,n_jobs=jobs, min_samples_split=1,class_weight="balanced",random_state=0,verbose=True)

    X=np.concatenate((clase0,clase1),0)
    Y=np.concatenate((label0,label1),0)

    clase0=None
    clase1=None
    label0=None
    label1=None

    print X.shape

    clf.fit(X,Y)
    if not exists(join(subject.dir,'classifiers_csfext')):
        makedirs(join(subject.dir,'classifiers_csfext'))

    if not exists(join(subject.dir,'classifiers_csfext',subject.name)):
        makedirs(join(subject.dir,'classifiers_csfext',subject.name))

    if not exists(join(subject.dir,'classifiers_csfext',subject.name,'classifier_'+str(n_estim)+'_'+str(depth))):
        makedirs(join(subject.dir,'classifiers_csfext',subject.name,'classifier_'+str(n_estim)+'_'+str(depth)))


    joblib.dump(clf, join(subject.dir,'classifiers_csfext',subject.name,'classifier_'+str(n_estim)+'_'+str(depth),'clasifier_csf_ext_'+str(n_estim)+'_'+str(depth)+'.pkl'))

    subject.classifier_ext_csf = 1
    subject.classifier_ext_csf_path = join(subject.dir,'classifiers_csfext',subject.name,'classifier_'+str(n_estim)+'_'+str(depth),'clasifier_csf_ext_'+str(n_estim)+'_'+str(depth)+'.pkl')

    return subject

def train_csfext(dir,subject_list,n_estim,depth,jobs):

    lista_0=[]
    lista_1=[]
    for subject_to_train in subject_list:


        print subject_to_train.name


        lista_0.append(subsampling(subject_to_train.class0_extcsf_descriptor_path,0.1))
        lista_1.append(subsampling(subject_to_train.class1_extcsf_descriptor_path, 0.3))



    print 'concatenando'
    clase0=np.concatenate(lista_0,0)
    lista_0=None


    clase1=np.concatenate(lista_1,0)

    lista_1=None

    print clase0.shape
    print clase1.shape

    label0=np.ones([clase0.shape[0],],dtype=int)

    label1=2*np.ones((clase1.shape[0],),dtype=int)




    print 'subset creado comenzando entrenamiento'


    clf = RandomForestClassifier(n_estimators=n_estim, max_depth=depth,n_jobs=jobs, min_samples_split=1,class_weight="balanced",random_state=0,verbose=True)

    X=np.concatenate((clase0,clase1),0)
    Y=np.concatenate((label0,label1),0)

    clase1=None
    clase2=None
    label1=None
    label2=None

    print X.shape

    clf.fit(X,Y)
    if not exists(join(dir,'classifiers_server_csfext')):
        makedirs(join(dir,'classifiers_server_csfext'))


    if not exists(join(dir,'classifiers_server_csfext','classifier_'+str(n_estim)+'_'+str(depth))):
        makedirs(join(dir,'classifiers_server_csfext','classifier_'+str(n_estim)+'_'+str(depth)))


    joblib.dump(clf, join(dir,'classifiers_server_csfext','classifier_'+str(n_estim)+'_'+str(depth),'clasifier_csf_ext_'+str(n_estim)+'_'+str(depth)+'.pkl'))

def test_csfext(subject,n_estim, depth,flag):



    print 'creating csf ext'
    if flag == 0:

        clase1 = np.load(subject.clase0_extcsf_descriptor_path, 'r')
        clase2 = np.load(subject.clase1_extcsf_descriptor_path, 'r')
        indices_1 = np.load(subject.clase0_extcsf_indices_path, 'r')
        indices_2 = np.load(subject.clase1_extcsf_indices_path, 'r')
        indices = np.concatenate((indices_1, indices_2), 0)
        data_test = np.concatenate((clase1, clase2), 0)

    if flag == 1:

        data_test = np.load(subject.test_extcsf_descriptor_path, 'r')
        indices = np.load(subject.test_extcsf_indices_path, 'r')

    if flag == 2:

        data_test = np.load(subject.test_extcsf_descriptor_path, 'r')
        indices = np.load(subject.test_extcsf_indices_path, 'r')



    if flag > 2:
        'error invalid flag'
        return -1


    mask_for_metadata,metadata_header = load(subject.brain_mask_path)

    result = np.zeros(mask_for_metadata.shape, dtype=float)


    print 'Concatenate done'

    clase1 = None
    clase2 = None

    test_split = np.array_split(data_test, 10, 0)
    if flag != 2:
        clf = joblib.load(join(subject.dir,'classifiers_csfext',subject.name,'classifier_'+str(n_estim)+'_'+str(depth),'clasifier_csf_ext_'+str(n_estim)+'_'+str(depth)+'.pkl'))
    else:
        clf = joblib.load(subject.classifier_ext_csf_path)

    pred_labels_all = np.array([], dtype=np.int).reshape(0, )

    for j in test_split:
        print  j.shape

        pred_labels = clf.predict_proba(j)


        pred_labels_all = np.concatenate((pred_labels_all, pred_labels[:, 1]), 0)

    contador = 0

    for pred_value in pred_labels_all:
        result[indices[contador][0]][indices[contador][1]][indices[contador][2]] = pred_value

        contador = contador + 1


    result = (result > 0.5)
    if flag == 0:
        print 'saving...'

        if not exists(join(subject.dir, 'results_csfext')):
            makedirs(join(subject.dir, 'results_csfext'))

        if not exists(join(subject.dir, 'results_csfext', subject.name)):
            makedirs(join(subject.dir, 'results_csfext', subject.name))

        if not exists(
                join(subject.dir, 'results_csfext', subject.name, 'results_csfext_' + str(n_estim) + '_' + str(depth))):
            makedirs(
                join(subject.dir, 'results_csfext', subject.name, 'results_csfext_' + str(n_estim) + '_' + str(depth)))


        save(result, join(subject.dir, 'results_csfext', subject.name, 'results_csfext_' + str(n_estim) + '_' + str(depth), 'csf_ext.nii.gz'),metadata_header)
        subject.csf_ext= 1
        subject.csf_ext_path = join(subject.dir, 'results_csfext', subject.name, 'results_csfext_' + str(n_estim) + '_' + str(depth), 'csf_ext.nii.gz')

        print 'saved'

    if flag == 1:
        print 'saving...'

        if not exists(join(subject.dir, 'results_csfext')):
            makedirs(join(subject.dir, 'results_csfext'))

        if not exists(join(subject.dir, 'results_csfext', subject.name)):
            makedirs(join(subject.dir, 'results_csfext', subject.name))

        if not exists(
                join(subject.dir, 'results_csfext', subject.name, 'test_results_csfext_' + str(n_estim) + '_' + str(depth))):
            makedirs(
                join(subject.dir, 'results_csfext', subject.name, 'test_results_csfext_' + str(n_estim) + '_' + str(depth)))

        save(result, join(subject.dir, 'results_csfext', subject.name, 'test_results_csfext_' + str(n_estim) + '_' + str(depth),
                          'csf_ext.nii.gz'), metadata_header)

        subject.csf_ext= 1
        subject.csf_ext_path = join(subject.dir, 'results_csfext', subject.name, 'test_results_csfext_' + str(n_estim) + '_' + str(depth), 'csf_ext.nii.gz')

        print 'saved'
    if flag == 2:
        print 'saving...'

        if not exists(join(subject.intermediate_path, 'results_csf_ext')):
            makedirs(join(subject.intermediate_path, 'results_csf_ext'))


        if not exists(
                join(subject.intermediate_path, 'results_csf_ext', 'results_csfext_' + str(n_estim) + '_' + str(depth))):
            makedirs(
                join(subject.intermediate_path, 'results_csf_ext', 'results_csfext_' + str(n_estim) + '_' + str(depth)))

        save(result, join(subject.intermediate_path, 'results_csf_ext', 'results_csfext_' + str(n_estim) + '_' + str(depth),
                          'csf_ext.nii.gz'), metadata_header)

        subject.csf_ext= 1
        subject.csf_ext_path = join(subject.intermediate_path, 'results_csf_ext','results_csfext_' + str(n_estim) + '_' + str(depth), 'csf_ext.nii.gz')

        print 'saved'
    return subject

def features_importance_csfext(subject,n_estim,depth):



    clf = joblib.load(join(subject.dir,'classifiers_csfext',subject.name,'classifier_'+str(n_estim)+'_'+str(depth),'clasifier_csf_ext_'+str(n_estim)+'_'+str(depth)+'.pkl'))


    print clf.feature_importances_ * 100

    plt.plot(clf.feature_importances_)



    plt.show()

def create_features_ms(subject,flag):


    print "------------Creating MS features-----------------"

    mask, mask_header = load(subject.brain_mask_path)
    pv1, pv1_header = load(subject.GM_path)
    pv2, pv2_header = load(subject.WM_path)
    pv0, pv0_header = load(subject.CSF_path)
    flair, fl_header = load(subject.FLAIR_irs_path)
    mask_threshold=gray_matter_threshold(flair,(pv1>pv2)*(pv1>pv0))


    if flag == 0:
        mask_class, mask_class_header = load(subject.gt_ms_path)
        data_name = 'lesion'
    if flag == 1:
        mask_class, mask_class_header = load(subject.gt_ms_path)
        mask_class = (mask_class == 0)
        data_name = 'sane'


    if flag == 2:
        mask_class = 1
        data_name = 'brain'

    if flag > 2:

        'error invalid flag'
        return -1

    mask_voxels = mask * mask_threshold * mask_class

    iter = 0

    for vol_path in subject.return_irs_sequences():


        print 'Generating tissue and distance features---< ' +  '  ' + vol_path

        image_data, image_header = load(vol_path)

        spacing = medpy.io.header.get_pixel_spacing(image_header)

        intensities = medpy.features.intensity.intensities(image_data, mask_voxels)

        gaussian1 = medpy.features.intensity.local_mean_gauss(image_data, 3, spacing, mask_voxels)
        gaussian2 = medpy.features.intensity.local_mean_gauss(image_data, 5, spacing, mask_voxels)
        gaussian3 = medpy.features.intensity.local_mean_gauss(image_data, 7, spacing, mask_voxels)

        kernel1 = kernel_subs_creator(3, spacing)
        subs1 = medpy.features.intensity.intensities(ndimage.filters.convolve(image_data, kernel1), mask_voxels)
        kernel2 = kernel_subs_creator(5, spacing)
        subs2 = medpy.features.intensity.intensities(ndimage.filters.convolve(image_data, kernel2), mask_voxels)
        kernel3 = kernel_subs_creator(7, spacing)
        subs3 = medpy.features.intensity.intensities(ndimage.filters.convolve(image_data, kernel3), mask_voxels)

        if iter == 0:
            joined = medpy.features.utilities.join(intensities, gaussian1, gaussian2, gaussian3, subs1, subs2, subs3)
        else:
            joined = medpy.features.utilities.join(joined, intensities, gaussian1, gaussian2, gaussian3, subs1, subs2,
                                                   subs3)
        iter += 1
    print 'Generating tissue and distance features'
    spacing_indices = (1, 1, 1)

    spacing = medpy.io.header.get_pixel_spacing(fl_header)

    indices = medpy.features.indices(flair, spacing_indices, mask_voxels)

    distances_0 = medpy.features.intensity.centerdistance_xdminus1(flair, 0, spacing, mask_voxels)
    distances_1 = medpy.features.intensity.centerdistance_xdminus1(flair, 1, spacing, mask_voxels)
    distances_2 = medpy.features.intensity.centerdistance_xdminus1(flair, 2, spacing, mask_voxels)




    spacing = medpy.io.header.get_pixel_spacing(fl_header)


    ext_csf, ext_csf_header = load(subject.csf_ext_path)
    #
    dist_transform = ndimage.distance_transform_edt(ext_csf == 0, sampling=spacing)
    dist_transform_feature = medpy.features.intensity.intensities(dist_transform, mask_voxels)


    intensities_pv0 = medpy.features.intensity.intensities(pv0, mask_voxels)
    intensities_pv1 = medpy.features.intensity.intensities(pv1, mask_voxels)
    intensities_pv2 = medpy.features.intensity.intensities(pv2, mask_voxels)

    gaussian1_pv0 = medpy.features.intensity.local_mean_gauss(pv0, 3, spacing, mask_voxels)
    gaussian2_pv0 = medpy.features.intensity.local_mean_gauss(pv0, 7, spacing, mask_voxels)
    gaussian3_pv0 = medpy.features.intensity.local_mean_gauss(pv0, 15, spacing, mask_voxels)

    gaussian1_pv1 = medpy.features.intensity.local_mean_gauss(pv1, 3, spacing, mask_voxels)
    gaussian2_pv1 = medpy.features.intensity.local_mean_gauss(pv1, 7, spacing, mask_voxels)
    gaussian3_pv1 = medpy.features.intensity.local_mean_gauss(pv1, 15, spacing, mask_voxels)

    gaussian1_pv2 = medpy.features.intensity.local_mean_gauss(pv2, 3, spacing, mask_voxels)
    gaussian2_pv2 = medpy.features.intensity.local_mean_gauss(pv2, 7, spacing, mask_voxels)
    gaussian3_pv2 = medpy.features.intensity.local_mean_gauss(pv2, 15, spacing, mask_voxels)

    joined = medpy.features.utilities.join(joined, distances_0, distances_1, distances_2, intensities_pv0, intensities_pv1,
                                           intensities_pv2, gaussian1_pv0, gaussian2_pv0, gaussian3_pv0, gaussian1_pv1
                                           , gaussian2_pv1, gaussian3_pv1, gaussian1_pv2, gaussian2_pv2, gaussian3_pv2,
                                            dist_transform_feature)












    if flag == 0:
        if not exists(join(subject.dir, 'descriptors')):
            makedirs(join(subject.dir, 'descriptors'))

        np.save(join(subject.dir, 'descriptors', data_name + '_descriptor'), joined)

        np.save(join(subject.dir, 'descriptors', data_name + '_indices'), indices)
        subject.class0_ms_descriptor_path = join(subject.dir, 'descriptors', data_name + '_descriptor.npy')
        subject.class0_ms_indices_path = join(subject.dir, 'descriptors', data_name + '_indices.npy')
    if flag == 1:
        if not exists(join(subject.dir, 'descriptors')):
            makedirs(join(subject.dir, 'descriptors'))

        np.save(join(subject.dir, 'descriptors', data_name + '_descriptor'), joined)

        np.save(join(subject.dir, 'descriptors', data_name + '_indices'), indices)
        subject.class1_ms_descriptor_path = join(subject.dir, 'descriptors', data_name + '_descriptor.npy')
        subject.class1_ms_indices_path = join(subject.dir, 'descriptors',  data_name + '_indices.npy')

    if flag == 2:
        if not exists(join(subject.intermediate_path, 'descriptors')):
            makedirs(join(subject.intermediate_path, 'descriptors'))

        np.save(join(subject.intermediate_path, 'descriptors', data_name + '_descriptor'), joined)

        np.save(join(subject.intermediate_path, 'descriptors', data_name + '_indices'), indices)
        subject.test_ms_descriptor_path = join(subject.intermediate_path, 'descriptors',data_name + '_descriptor.npy')
        subject.test_ms_indices_path = join(subject.intermediate_path, 'descriptors', data_name + '_indices.npy')
    print 'ms features saved'
    return subject

def train_subject_ms(subject,subject_list,n_estim,depth,jobs):

    lista_0=[]
    lista_1=[]
    for subject_to_train in subject_list:

        if subject_to_train.name!=subject.name:
            print subject_to_train.name


            lista_0.append(subsampling(subject_to_train.class0_ms_descriptor_path,1))
            lista_1.append(subsampling(subject_to_train.class1_ms_descriptor_path, 1))




    print 'concatenando'
    clase0=np.concatenate(lista_0,0)
    lista_0=None


    clase1=np.concatenate(lista_1,0)

    lista_1=None

    print clase0.shape
    print clase1.shape

    label0=2*np.ones([clase0.shape[0],],dtype=int)

    label1=np.ones((clase1.shape[0],),dtype=int)




    print 'subset creado comenzando entrenamiento'


    clf = RandomForestClassifier(n_estimators=n_estim, max_depth=depth,n_jobs=jobs, min_samples_split=1,class_weight="balanced_subsample",random_state=0,verbose=True)

    X=np.concatenate((clase0,clase1),0)
    Y=np.concatenate((label0,label1),0)

    clase0=None
    clase1=None
    label0=None
    label1=None

    print X.shape

    clf.fit(X,Y)
    if not exists(join(subject.dir,'classifiers_ms')):
        makedirs(join(subject.dir,'classifiers_ms'))

    if not exists(join(subject.dir,'classifiers_ms',subject.name)):
        makedirs(join(subject.dir,'classifiers_ms',subject.name))

    if not exists(join(subject.dir,'classifiers_ms',subject.name,'classifier_'+str(n_estim)+'_'+str(depth))):
        makedirs(join(subject.dir,'classifiers_ms',subject.name,'classifier_'+str(n_estim)+'_'+str(depth)))


    joblib.dump(clf, join(subject.dir,'classifiers_ms',subject.name,'classifier_'+str(n_estim)+'_'+str(depth),'clasifier_ms_'+str(n_estim)+'_'+str(depth)+'.pkl'))


    subject.classifier_ms = 1
    subject.classifier_ms_path = join(subject.dir,'classifiers_ms',subject.name,'classifier_'+str(n_estim)+'_'+str(depth),'clasifier_ms_'+str(n_estim)+'_'+str(depth)+'.pkl')

    return subject

def train_ms(dir,subject_list,n_estim,depth,jobs):

    lista_0=[]
    lista_1=[]
    for subject_to_train in subject_list:


        print subject_to_train.name


        lista_0.append(subsampling(subject_to_train.class0_ms_descriptor_path,1))
        lista_1.append(subsampling(subject_to_train.class1_ms_descriptor_path, 1))



    print 'concatenando'
    clase1=np.concatenate(lista_0,0)
    lista_0=None


    clase2=np.concatenate(lista_1,0)

    lista_1=None

    print clase1.shape
    print clase2.shape

    label1=2*np.ones([clase1.shape[0],],dtype=int)

    label2=np.ones((clase2.shape[0],),dtype=int)




    print 'subset creado comenzando entrenamiento'


    clf = RandomForestClassifier(n_estimators=n_estim, max_depth=depth,n_jobs=jobs, min_samples_split=1,class_weight="balanced",random_state=0,verbose=True)

    X=np.concatenate((clase1,clase2),0)
    Y=np.concatenate((label1,label2),0)

    clase1=None
    clase2=None
    label1=None
    label2=None

    print X.shape

    clf.fit(X,Y)
    if not exists(join(dir, 'classifiers_server_ms')):
        makedirs(join(dir, 'classifiers_server_ms'))

    if not exists(join(dir, 'classifiers_server_ms', 'classifier_' + str(n_estim) + '_' + str(depth))):
        makedirs(join(dir, 'classifiers_server_ms', 'classifier_' + str(n_estim) + '_' + str(depth)))

    joblib.dump(clf, join(dir, 'classifiers_server_ms', 'classifier_' + str(n_estim) + '_' + str(depth),
                          'clasifier_ms_' + str(n_estim) + '_' + str(depth) + '.pkl'))

def test_ms(subject,n_estim, depth,flag):



    print 'creating ms prob mask'
    if flag == 0:

        clase1 = np.load(subject.class0_ms_descriptor_path, 'r')
        clase2 = np.load(subject.class1_ms_descriptor_path, 'r')
        indices_1 = np.load(subject.class0_ms_indices_path, 'r')
        indices_2 = np.load(subject.class1_ms_indices_path, 'r')
        indices = np.concatenate((indices_1, indices_2), 0)
        data_test = np.concatenate((clase1, clase2), 0)

    if flag == 1:

        data_test = np.load(subject.test_ms_descriptor_path, 'r')
        indices = np.load(subject.test_ms_indices_path, 'r')

    if flag == 2:
            data_test = np.load(subject.test_ms_descriptor_path, 'r')
            indices = np.load(subject.test_ms_indices_path, 'r')

    mask_for_metadata, metadata_header = load(subject.brain_mask_path)

    result = np.zeros(mask_for_metadata.shape, dtype=float)


    print 'Concatenate done'

    clase1 = None
    clase2 = None

    test_split = np.array_split(data_test, 10, 0)

    if flag != 2:
        clf = joblib.load(
            join(subject.dir, 'classifiers_ms', subject.name, 'classifier_' + str(n_estim) + '_' + str(depth),
                 'clasifier_ms_' + str(n_estim) + '_' + str(depth) + '.pkl'))
    else:
        clf = joblib.load(subject.classifier_ms_path)

    pred_labels_all = np.array([], dtype=np.int).reshape(0, )

    for j in test_split:
        print  j.shape

        pred_labels = clf.predict_proba(j)


        pred_labels_all = np.concatenate((pred_labels_all, pred_labels[:, 1]), 0)

    contador = 0

    for pred_value in pred_labels_all:
        result[indices[contador][0]][indices[contador][1]][indices[contador][2]] = pred_value

        contador = contador + 1


    if flag == 0:
        print 'saving...'

        if not exists(join(subject.dir, 'results_ms')):
            makedirs(join(subject.dir, 'results_ms'))

        if not exists(join(subject.dir, 'results_ms', subject.name)):
            makedirs(join(subject.dir, 'results_ms', subject.name))

        if not exists(
                join(subject.dir, 'results_ms', subject.name, 'results_ms_' + str(n_estim) + '_' + str(depth))):
            makedirs(
                join(subject.dir, 'results_ms', subject.name, 'results_ms_' + str(n_estim) + '_' + str(depth)))


        save(result, join(subject.dir, 'results_ms', subject.name, 'results_ms_' + str(n_estim) + '_' + str(depth), 'lesion_prob.nii.gz'),metadata_header)
        subject.lesion_prob = 1
        subject.lesion_prob_path = join(subject.dir, 'results_ms', subject.name, 'results_ms_' + str(n_estim) + '_' + str(depth), 'lesion_prob.nii.gz')
        print 'saved'

    if flag == 1:
        print 'saving...'

        if not exists(join(subject.dir, 'results_ms')):
            makedirs(join(subject.dir, 'results_ms'))

        if not exists(join(subject.dir, 'results_ms', subject.name)):
            makedirs(join(subject.dir, 'results_ms', subject.name))

        if not exists(
                join(subject.dir, 'results_ms', subject.name, 'test_results_ms_' + str(n_estim) + '_' + str(depth))):
            makedirs(
                join(subject.dir, 'results_ms', subject.name, 'test_results_ms_' + str(n_estim) + '_' + str(depth)))

        save(result, join(subject.dir, 'results_ms', subject.name, 'test_results_ms_' + str(n_estim) + '_' + str(depth),
                          'lesion_prob.nii.gz'), metadata_header)
        subject.lesion_prob = 1
        subject.lesion_prob_path = join(subject.dir, 'results_ms', subject.name,
                                        'test_results_ms_' + str(n_estim) + '_' + str(depth),
                                        'lesion_prob.nii.gz')

        print 'saved'

    if flag == 2:
        print 'saving...'

        if not exists(join(subject.intermediate_path, 'results_rf')):
            makedirs(join(subject.intermediate_path, 'results_rf'))



        if not exists(
                join(subject.intermediate_path, 'results_rf',  'test_results_rf_' + str(n_estim) + '_' + str(depth))):
            makedirs(
                join(subject.intermediate_path, 'results_rf',  'test_results_rf_' + str(n_estim) + '_' + str(depth)))

        save(result, join(subject.intermediate_path, 'results_rf', 'test_results_rf_' + str(n_estim) + '_' + str(depth),
                          'lesion_prob.nii.gz'), metadata_header)
        subject.lesion_prob = 1
        subject.lesion_prob_path = join(subject.intermediate_path, 'results_rf',
                                        'test_results_rf_' + str(n_estim) + '_' + str(depth),
                                        'lesion_prob.nii.gz')

        print 'saved'


    return subject

def lesion_growing(subject,theta,beta_grow,flag=0):
    print "-----------------Creating final lesion mask----------------------"
    print('Grow start')
    mask, mask_header = load(subject.brain_mask_path)
    pv1, pv1_header = load(subject.GM_path)
    pv2, pv2_header = load(subject.WM_path)
    pv0, pv0_header = load(subject.CSF_path)
    flair, fl_header = load(subject.FLAIR_irs_path)

    mask_for_metadata, metadata_header = load(subject.brain_mask_path)

    lesiones_prob,lesiones_header_prob = load(subject.lesion_prob_path)

    mask_threshold=gray_matter_threshold(flair,(pv1>pv2)*(pv1>pv0))

    lesiones_clas = (lesiones_prob > theta)*mask_threshold*mask
    lesiones_clas=filter.binary.size_threshold(lesiones_clas,10, comp='lt', structure=None)

    pv0_corrected = (pv0>pv1) * (pv0>pv2) * (lesiones_clas==0)*mask
    pv1_corrected = (pv1>pv2) * (pv1>pv0) * (lesiones_clas==0)* mask
    pv2_corrected = (pv2>pv1) * (pv2>pv0) * (lesiones_clas==0)*mask


    csf = medpy.features.intensities(flair,pv0_corrected)
    GM = medpy.features.intensities(flair,pv1_corrected)
    WM = medpy.features.intensities(flair,pv2_corrected)
    les = medpy.features.intensities(flair,lesiones_clas)

    brain_data = np.concatenate((csf,GM,WM))



    gaussMixture= mixture.GMM(n_components=3)
    gaussMixture.means_ = np.asarray([[np.mean(csf)],[np.mean(GM)],[np.mean(WM)]])
    gmm_data = np.array([brain_data]).transpose()



    gaussMixture.fit(gmm_data)

    if len(les)>10:
        alpha, loc, beta=stats.gamma.fit(les)


    struct2 = ndimage.generate_binary_structure(3, 1)

    grow_counter =100000000
    while grow_counter >300:

        spacing_indices=(1,1,1)

        les_grow = binary_dilation(lesiones_clas,structure= struct2, iterations=1)



        contador=0

        indices=medpy.features.indices(flair,spacing_indices, les_grow*(lesiones_clas==0))
        resultado=np.zeros(lesiones_clas.shape,dtype=float)
        kernel = np.ones((3,3,3))/(3*3*3)

        N_les=ndimage.filters.convolve(lesiones_prob,kernel)
        posible_les=medpy.features.intensity.intensities( N_les,  les_grow*(lesiones_clas==0))

        for pred_value in posible_les:

            flair_valor = flair[indices[contador][0]][indices[contador][1]][indices[contador][2]]

            resultado[indices[contador][0]][indices[contador][1]][indices[contador][2]]=(stats.gamma.pdf(flair_valor,alpha, loc=loc, scale=beta)*pred_value)/(np.exp(gaussMixture.score_samples(flair_valor)[0])*(1-pred_value))

            contador=contador+1

        lesiones_clas = (resultado>beta_grow) +lesiones_clas
        grow_counter = np.sum(resultado>beta_grow)


    print('Grow done')
    if flag==0:
        if not exists(join(subject.dir, 'results_grow')):
            makedirs(join(subject.dir, 'results_grow'))


        if not exists(
                join(subject.dir, 'results_grow', 'results_grow_' + str(int(theta*100)) + '_' + str(int(beta_grow*100)) )):
            makedirs(
                join(subject.dir, 'results_grow',  'results_grow_' + str(int(theta*100)) + '_' + str(int(beta_grow*100))))

        save(lesiones_clas,join(subject.dir, 'results_grow','results_grow_' + str(int(theta*100)) + '_' + str(int(beta_grow*100)), 'lesion_mask.nii.gz'),metadata_header)

        subject.lesion_mask_final = 1
        subject.lesion_mask_final_path = join(subject.dir, 'results_grow', 'results_grow_' + str(int(theta*100)) + '_' + str(int(beta_grow*100)),'lesion_mask.nii.gz')

    if flag==1:



        if not exists(
                join(subject.dir, 'results' )):
            makedirs(
                join(subject.dir, 'results'))

        save(lesiones_clas,join(subject.dir, 'results', 'lesion_mask.nii.gz'),metadata_header)

        subject.lesion_mask_final = 1
        subject.lesion_mask_final_path = join(subject.dir, 'results','lesion_mask.nii.gz')
        print 'Lesion mask generated  '


    return subject

def calc_result(subject_list, xls_name):

    workbook = xlsxwriter.Workbook(xls_name)


    worksheet = workbook.add_worksheet()

    # Some data we want to write to the worksheet.
    titulo = (
    'Dice ',
    'Precision',
    'FPR',
    'TPR',
    'ASSD',
    'Densidad'
    )

    # Start from the first cell. Rows and columns are zero indexed.
    row = 3
    col = 3
    col_inicio=col

    # Iterate over the data and write it out row by row.
    for item in (titulo):
        worksheet.write(row, col,     item)

        col += 1

    for subject in subject_list:

        row+=1
        col=col_inicio-1


        lesiones, lesiones_header = load(subject.gt_ms_path)
        mask, mask_header = load(subject.brain_mask_path)

        lesiones_clas,lesiones_header_clas = load(subject.lesion_prob_path)
        lesiones_clas = lesiones_clas > 0.5

        flair, fl_header = load(subject.FLAIR_irs_path)
        spacing = medpy.io.header.get_pixel_spacing(fl_header)
        tpr1 = float(np.sum(lesiones_clas*lesiones))/np.sum(lesiones)
        fpr1 = np.sum((np.ones(lesiones.shape)-lesiones)*(lesiones_clas))/np.sum(lesiones==1)

        dc= medpy.metric.dc(lesiones,lesiones_clas)
        recall=medpy.metric.recall(lesiones_clas,lesiones)
        precision=medpy.metric.precision(lesiones_clas,lesiones)
        densidad = np.sum(lesiones_clas)/np.sum(mask)
        assd = medpy.metric.assd(lesiones,lesiones_clas,voxelspacing=spacing)
        nombre= subject.name

        print  nombre
        print "dc = " + str(dc)
        print "recall = " + str(recall)
        print "precision = " + str(precision)
        print "TPR_voxel = "+ str(tpr1)
        print "FPS_voxel = "+ str(fpr1)
        print "assd_voxel = "+ str(assd)
        print "densidad = "+ str(densidad)


        datos = (
        nombre,
        dc,
        precision,
        fpr1,
        tpr1,
        assd,
        densidad
        )

        for item in (datos):
            worksheet.write(row, col,item)


            col += 1






    col=col_inicio-1
    worksheet.write(row, col,'medias')
    worksheet.write_formula(row, col+1,'=AVERAGE(D$5:D$19)')
    worksheet.write_formula(row, col+2,'=AVERAGE(E$5:E$19)')
    worksheet.write_formula(row, col+3,'=AVERAGE(F$5:F$19)')
    worksheet.write_formula(row, col+4,'=AVERAGE(G$5:G$19)')
    worksheet.write_formula(row, col+5,'=AVERAGE(H$5:H$19)')
    worksheet.write_formula(row, col+5,'=AVERAGE(I$5:I$19)')
    workbook.close()

def gunzip_T1_pp(subject):

    vol, header = load(subject.T1_pp_path)
    save(vol, join(subject.intermediate_path,'T1_pp.nii'), header)
    subject.T1_pp_gunzip = 1
    subject.T1_pp_gunzip_path=join(subject.intermediate_path, 'T1_pp.nii')
    return subject


def gunzip_T1_unpp(subject):


    vol, header = load(subject.T1_path)
    save(vol, join(subject.intermediate_path, 'T1.nii'), header)
    subject.T1_gunzip = 1
    subject.T1_gunzip_path = join(subject.intermediate_path, 'T1.nii')
    return subject