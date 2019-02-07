import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm 
from matplotlib import mlab
from scipy.ndimage.filters import gaussian_filter
import argparse

def get_args():
    parser = argparse.ArgumentParser(
    "Program to Generate Petroleomic Data from Formula ID results output by Bruker DataAnalysis SmartFormula Routine")
    parser.add_argument('-data', help = 'Path to .csv data file with Formula ID results')
    parser.add_argument('--err_plots', help = "1 = Generate scatter and histogram plots, 0 = skip generation", 
                        default = '0')
    parser.add_argument('--w', help = 'Write non-plot data to .txt file, 1 == True, 0 == False', default = '0')
    parser.add_argument('--dbe_plots', help = "1 = Generate DBE vs C number plots, 0 = skip generation", 
                        default = '0')
    parser.add_argument('--class_plots', help = "1 = Generate class distribution plots for a standard set of classes, 0 = skip generation", 
                        default = '0')
    parser.add_argument('--keep_class', help =
     'If provided, this generates additional class dist plot based on a' + 
     'list of comma separated classes to keep EX. --keep_class NO2,O2,N3', default = None)
    parser.add_argument('--comp_data', help = 'Path to addtional dataset. ' + 
    'If provided, generates a class plot comparing two compound identificiation results', default = None)
    parser.add_argument('--comp_keep_class', help = 'Classes to keep for class distribution comparison\
         entered as comma separated EX. --comp_keep_class NO2,O2,S', default = None)
    

    args = parser.parse_args()


    return args

def get_carbon_num(formula):
    if 'H' in formula:
        s = formula.split('H')[0]
        if len(s) == 1:
            num = 1
        else:
            num = int(s.split('C')[1])
    elif 'N' in formula:
        s = formula.split('N')[0]
        if len(s) == 1:
            num = 1
        else:
            num = int(s.split('C')[1])
    elif 'O' in formula:
        s = formula.split('O')[0]
        if len(s) == 1:
            num = 1
        else:
            num = int(s.split('C')[1])
    elif 'S' in formula:
        s = formula.split('S')[0]
        if len(s) == 1:
            num = 1
        else:
            num = int(s.split('C')[1])
    return num


def extract_c_nums(df):
    l = df.tolist()
    carbon_nums = []
    for e in l:
        num = get_carbon_num(e)
        carbon_nums.append(num)
    return carbon_nums

def average_dbe(df):
    I = df['I']/df['I'].max()
    avg_dbe = np.sum(np.array(df['rdb'])*I)/np.sum(I)
    return avg_dbe


def plot_ppm_error(measured, reference, file, bins = 10, save = False):
    errors = ((measured - reference)/reference)*1e6
    plt.figure()
    plt.hist(errors, bins = bins, range = (-1,1), color = 'r')
    plt.tick_params(
        axis='y',         
        which='both',      # both major and minor ticks are affected
        left=False,      # ticks along the bottom edge are off
        right=False,         # ticks along the top edge are off
        labelleft=False) # labels along the bottom edge are off
    plt.tick_params(axis = 'x', labelsize =15)
    plt.xlabel("ppm mass error", fontsize = 'x-large')
    plt.gcf().subplots_adjust(bottom=0.15)
    if save == True:
        save_path = 'hist_' + file.split('.')[0] + '.png'
        plt.savefig(save_path, dpi = 1200)
        print("Histogram saved as " + save_path)
    rms_err = np.sqrt((np.sum(errors**2))/len(errors))
    sdev_err = np.std(errors)
    mean_err = np.mean(errors)
    nc = len(errors)
    return errors, nc, mean_err, sdev_err, rms_err

def plot_scatter_ppm(measured, errors, file):
    plt.figure()
    plt.scatter(measured, errors, s = 5, c = 'k')
    plt.tick_params(axis = 'x', labelsize =15)
    plt.xlabel("Measured m/z", fontsize = 'xx-large')
    plt.tick_params(axis = 'y', labelsize =15)
    plt.ylabel("ppm mass error", fontsize = 'xx-large')
    plt.ylim((-1,1))
    plt.xlim((np.min(measured) - 10.0, np.max(measured) + 10.0))
    plt.gcf().subplots_adjust(bottom=0.15, left = 0.15 )
    save_path = 'scatter_' + file.split('.')[0] + '.png'
    plt.savefig(save_path, dpi = 1200)
    print("Scatter plot saved as " + save_path)

def plot_DBEvC(df, carbon_nums, file):
    dbe = np.array(df['rdb'])
    carbon_nums = np.array(carbon_nums)
    fig = plt.figure()
    fig.set_size_inches(5.0, 4.0)
    plt.hexbin(carbon_nums, dbe, cmap=cm.jet, gridsize = (50,50))
    plt.axis([carbon_nums.min(), carbon_nums.max(), dbe.min(), dbe.max()])
    plt.tick_params(axis = 'x', labelsize =15)
    plt.xlabel("Carbon Number", fontsize = 'x-large')
    plt.tick_params(axis = 'y', labelsize =15)
    plt.ylabel("DBE", fontsize = 'x-large')
    plt.gcf().subplots_adjust(left = 0.2, bottom=0.15)
    cb = plt.colorbar(shrink = 1, pad = 0.025)
    cb.set_label('Occurrence', fontsize = 'x-large')
    dbe_path = 'DBEvsC_' + file.split('.')[0] + '.png'
    plt.savefig(dbe_path, dpi = 1200)
    print("DBE vs C plot saved as " + dbe_path)
    
def new_ticks(given_locs, edges):
    # Returns new ticks and labels for a matplotlib tick input
    n = len(given_locs[1:])
    labels = list(np.linspace(0, edges[-1], n).round().astype('int').astype('str'))
    return given_locs[1:], labels

def smooth_hist2d(x, y, file, s=2, bins=[100,100]):
    heatmap, xedges, yedges = np.histogram2d(x, y, bins=bins)
    heatmap = gaussian_filter(heatmap, sigma=s)
    fig = plt.figure()
    fig.set_size_inches(5.0, 4.0)
    plt.imshow(heatmap.T, aspect = 'equal', origin='lower', cmap=cm.jet)
    xlocs, xlabels = new_ticks(plt.xticks()[0], xedges)
    ylocs, ylabels = new_ticks(plt.yticks()[0], yedges)
    plt.xticks(xlocs, xlabels)
    plt.yticks(ylocs, ylabels)
    plt.tick_params(axis = 'x', labelsize =15)
    plt.xlabel("Carbon Number", fontsize = 'x-large')
    plt.tick_params(axis = 'y', labelsize =15)
    plt.ylabel("DBE", fontsize = 'x-large')
    plt.gcf().subplots_adjust(bottom=0.15)
    cb = plt.colorbar(shrink = 0.9, pad = 0.025)
    cb.set_label('Rel. Abundance', fontsize = 'x-large')
    blur_path = 'DBEvsC_Blur' + file.split('.')[0] + '.png'
    #cb.set_label('Occurrence', fontsize = 'x-large')
    plt.savefig(blur_path, dpi = 1200)
    print('Gaussian blur saved as ' + blur_path)
    return None


def hc_dbe(df):
    df = df.reset_index()
    dbes = []
    ints = []
    for i in range(len(df)):
        formula = df['Ion Formula'][i]
        if ('N' in formula) or ('O' in formula) or ('S' in formula):
            continue
        else:
            dbes.append(df['rdb'][i])
            ints.append(df['I'][i])
    return np.array(dbes), np.array(ints)

def get_class(formula):
    if 'N' in formula:
        cl = 'N' + formula.split('N')[1]
    elif 'O' in formula:
        cl = 'O' + formula.split('O')[1]
    elif 'S' in formula:
        cl = 'S' + formula.split('S')[1]
    else:
        cl = 'HC'
    return cl

def assign_class(formulas):
    cl_list = []
    for f in formulas:
        cl = get_class(f)
        cl_list.append(cl)
    
    cl_list = np.array(cl_list)
    
    return cl_list    
    
def determine_class_abund(formulas, intensity, ref_list):
    cl_list = assign_class(formulas) # list of classes (np.array)
    unique_classes = np.unique(cl_list) # Unique assignments in sample
    # List corresponding to publishing standard by comparison with a reference list.
    pub_list = list(set(unique_classes) and set(ref_list)) 
    
    abund_bins = []
    for cl in pub_list:

        abund = np.sum(intensity*[cl_list == cl])
        abund_bins.append(abund)
    
    return np.array(abund_bins), np.array(pub_list)
    
def make_subscripts(labels):
    SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

    for i in range(len(labels)):
        if not labels[i].isalpha():
            labels[i] = labels[i].translate(SUB)
    return labels

def plot_class_dist(abund_bins, labels, file):
    width = 0.5
    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.set_size_inches(3.25, 2.45)
    rects1 = ax.bar(np.arange(len(abund_bins)), abund_bins, width, color='b')
    ax.set_ylabel('Relative Abundance', fontsize = 6)
    ax.set_xlabel('Class', fontsize = 6)
    ax.set_xticks(np.arange(len(labels)))
    plt.setp(ax.get_xticklabels(), rotation = 90, fontsize=4.5)
    plt.setp(ax.get_yticklabels(), fontsize=6)
    ax.set_xticklabels(list(labels))
    plt.gcf().subplots_adjust(bottom=0.30)
    plt.gcf().subplots_adjust(left=0.30)
    plt.savefig('classdist_' + file.split('.')[0] + '.png', dpi = 1200)

def reduce_classes(labels, abunds, keep_list):
    new_labels = []
    new_abunds = []
    for i in range(len(labels)):
        if labels[i] in keep_list:
            new_labels.append(labels[i])
            new_abunds.append(abunds[i])
    
    return new_labels, new_abunds

if __name__ == "__main__":


    args = get_args()
    file = args.data
    df = pd.read_csv(file)
    dfmzs = df[pd.notnull(df['Meas. m/z'])]
    measured = np.array(dfmzs['Meas. m/z'])
    reference = np.array(dfmzs['m/z'])
    intensity = np.array(dfmzs['I']/dfmzs['I'].sum())
    carbon_nums = np.array(extract_c_nums(dfmzs['Ion Formula']))
    avg_dbe = average_dbe(dfmzs)
    plot_errors = args.err_plots == '1'

    errors, num_compounds, mean_err, sdev_err, rms_err = plot_ppm_error(measured, 
                                                                        reference,
                                                                        file, 
                                                                        bins = 20, 
                                                                        save = plot_errors)
    if plot_errors:
        plot_scatter_ppm(measured, errors, file)

    
    dbes, ints = hc_dbe(dfmzs)
    I = ints/np.sum(ints)
    avg_hc_dbe = np.sum(dbes*I)/np.sum(I)
    
    lessthan1 = np.sum([np.abs(errors) <= 1.0])
    print("Mean error (ppm): ", mean_err)
    print("Standard dev (ppm): ", sdev_err)
    print("Number of compounds: ", num_compounds)
    print("RMS error(ppm): ", rms_err)
    print('Compounds with error < 1 ppm: ', lessthan1)
    print("Average HC Class DBE ", avg_hc_dbe)

    if args.w == '1':
        txt_file = 'results_' + file.split('.')[0] + '.txt'
        with open(txt_file, 'w') as f:
            f.write("Number of compounds: %d" % num_compounds)
            f.write('\n')
            f.write("Mean error (ppm): %.8f"  % mean_err)
            f.write('\n')
            f.write("Standard dev (ppm): %.5f" % sdev_err)
            f.write('\n')
            f.write("RMS error(ppm): %.5f" % rms_err)
            f.write('\n')
            f.write('Average DBE: %.4f' % avg_dbe)
            f.write('\n')
            f.write('Compounds with error < 1 ppm: %d' % lessthan1)
            f.write('\n')
            f.write("Average HC Class DBE %.3f" %avg_hc_dbe)
        print("Results saved as " + txt_file)

    if args.dbe_plots == '1':
        dbe = np.array(dfmzs['rdb'])
        plot_DBEvC(dfmzs, carbon_nums, file)
        smooth_hist2d(carbon_nums, dbe, file)

    # Taken from Smith et al. 2018 DOI: 10.1021/acs.analchem.7b04159
    class_list1 = ['HC', 'N', 'NO','NO2','NO3','NO4','NO5','NOS', 
              'NO2S', 'NO3S', 'NO4S', 'NOS2', 'NS', 'NS2',
             'N2', 'N2O', 'N2O2', 'N2S', 'N2OS', 'OS', 'O2S', 
              'O3S', 'O4S', 'O5S', 'OS2', 'O2S2', 'O3S2', 'O4S2', 
              'O', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7' ]
        
    formulas = np.array(dfmzs['Ion Formula'])
    abund_bins, labels = determine_class_abund(formulas, intensity, class_list1)
    inds = np.argsort(labels)[::1]
    abund_bins = abund_bins[inds]
    labels = labels[inds]
    labels = make_subscripts(labels)

    if args.class_plots == '1':
        plot_class_dist(abund_bins, labels, file)

    if args.keep_class is not None:
        keep_list = make_subscripts(args.keep_class.split(','))
        new_labels, new_abunds = reduce_classes(labels, abund_bins, keep_list)
        plot_class_dist(new_abunds, new_labels, 'reduced_' + file)
    
    if args.comp_keep_class is None and args.comp_data is None:
        print("No comparison class distribution plot was generated because a class keep list and datafile were not provided.")

    elif args.comp_data_keep is not None and args.comp_data is not None:
        file2 = args.comp_data
        df2 = pd.read_csv(file2)
        dfmzs2 = df2[pd.notnull(df2['Meas. m/z'])]

        intensity2 = np.array(dfmzs2['I']/dfmzs2['I'].sum())
        formulas2 = np.array(dfmzs2['Ion Formula'])
        abund_bins2, labels2 = determine_class_abund(formulas2, intensity2, class_list1)
        inds2 = np.argsort(labels2)[::1]
        abund_bins2 = abund_bins2[inds2]
        labels2 = labels2[inds2]
        labels2 = make_subscripts(labels2)
        keep_list = make_subscripts(args.comp_keep_class.split(','))
        labels, abund_bins = reduce_classes(labels, abund_bins, keep_list)
        labels2, abund_bins2 = reduce_classes(labels2, abund_bins2, keep_list)
        width = 0.3
        fig = plt.figure()
        ax = fig.add_subplot(111)
        fig.set_size_inches(3.25, 2.45)
        rects1 = ax.bar(np.arange(len(abund_bins)), abund_bins, width, color='b')
        rects2 = ax.bar(np.arange(len(abund_bins2)) + width, abund_bins2, width, color='r')
        ax.set_ylabel('Relative Abundance', fontsize = 6)
        ax.set_xlabel('Class', fontsize = 6)
        ax.set_xticks(np.arange(len(labels)) + width/2)
        plt.setp(ax.get_xticklabels(), rotation = 90, fontsize=4.5)
        plt.setp(ax.get_yticklabels(), fontsize=6)
        ax.set_xticklabels(list(labels))
        plt.gcf().subplots_adjust(bottom=0.30)
        plt.gcf().subplots_adjust(left=0.30)
        ax.legend((rects1[0], rects2[0]), ('Toluene', 'Anisole'))

        class_comp_path = 'classdist_COMP_' + file.split('.')[0] + '_' + file2.split('.')[0] + '.png'
        plt.savefig(class_comp_path, dpi = 1200)
        print("Comparison of class distribution figure saved as "+class_comp_path)