import glob
import pandas as pd
import numpy as np
from numpy import power, exp
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from scipy.optimize import curve_fit


def auto_aa(data,R1a,R1b,Kab,Kba):
    a_11 = R1a + Kab
    a_22 = R1b + Kba
    a_sum = a_11 + a_22
    a_diff = a_11 - a_22
    lamb_1 = 0.5*(a_sum + np.sqrt(power(a_diff, 2) + 4*Kab*Kba))
    lamb_2 = 0.5*(a_sum - np.sqrt(power(a_diff, 2) + 4*Kab*Kba))
    lamb_diff = lamb_1 - lamb_2

    return (((lamb_1 - a_11)*exp(-lamb_2*data)) - ((lamb_2 - a_11)*exp(-lamb_1*data)))/lamb_diff

def auto_bb(data,R1a,R1b,Kab,Kba):
    a_11 = R1a + Kab
    a_22 = R1b + Kba
    a_sum = a_11 + a_22
    a_diff = a_11 - a_22
    lamb_1 = 0.5*(a_sum + np.sqrt(power(a_diff, 2) + 4*Kab*Kba))
    lamb_2 = 0.5*(a_sum - np.sqrt(power(a_diff, 2) + 4*Kab*Kba))
    lamb_diff = lamb_1 - lamb_2

    return (((lamb_1 - a_22)*exp(-lamb_2*data)) - ((lamb_2 - a_22)*exp(-lamb_1*data)))/lamb_diff

def cross_ab(data,R1a,R1b,Kab,Kba):
    a_11 = R1a + Kab
    a_22 = R1b + Kba
    a_21 = -Kab
    a_sum = a_11 + a_22
    a_diff = a_11 - a_22
    lamb_1 = 0.5*(a_sum + np.sqrt(power(a_diff, 2) + 4*Kab*Kba))
    lamb_2 = 0.5*(a_sum - np.sqrt(power(a_diff, 2) + 4*Kab*Kba))
    lamb_diff = lamb_1 - lamb_2

    return ((a_21*exp(-lamb_1*data)) - (a_21*exp(-lamb_2*data)))/lamb_diff

def cross_ba(data,R1a,R1b,Kab,Kba):
    a_11 = R1a + Kab
    a_22 = R1b + Kba
    a_12 = -Kba
    a_sum = a_11 + a_22
    a_diff = a_11 - a_22
    lamb_1 = 0.5*(a_sum + np.sqrt(power(a_diff, 2) + 4*Kab*Kba))
    lamb_2 = 0.5*(a_sum - np.sqrt(power(a_diff, 2) + 4*Kab*Kba))
    lamb_diff = lamb_1 - lamb_2
    
    return ((a_12*exp(-lamb_1*data)) - (a_12*exp(-lamb_2*data)))/lamb_diff


res = ['2M', '10T', '36G', '54G', '56T', '68L']
time = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]

time_arr = np.array(time)
comboX = np.hstack((time_arr,time_arr,time_arr,time_arr))
no = len(time)
def comboFunc(comboData,R1a,R1b,Kab,Kba):
    d1 = comboData[:no]         # auto_aa data
    d2 = comboData[no:2*no]     # auto_bb data
    d3 = comboData[2*no:3*no]   # cross_ab data
    d4 = comboData[3*no:4*no]   # cross_ba data
    
    res1 = auto_aa(d1,R1a,R1b,Kab,Kba)
    res2 = auto_bb(d2,R1a,R1b,Kab,Kba)
    res3 = cross_ab(d3,R1a,R1b,Kab,Kba)
    res4 = cross_ba(d4,R1a,R1b,Kab,Kba)
    
    return np.hstack((res1, res2, res3, res4))


files = glob.glob("*_edited.txt")

# output data file
fo = open('final_results.txt', 'w+')
fo.write("RES\tR1a\tR1b\tkab\tkba\n")

for file in files:
    data_int = pd.read_csv(file, delim_whitespace=' ', index_col=0, usecols=[x for x in range(0,14)])
    print(file)
    # fo.write("{}\n".format(file))
    #print(data_int)
    auto_cross = []
    for a in res:
        auto_cross.append(a + 'a')
        auto_cross.append(a + 'b')
        auto_cross.append(a + 'ba')
        auto_cross.append(a + 'ab')

    data_dict = {}

    for p in auto_cross:
        ints = []
        for t in time:
            if t == 0 and (p[-2:] == 'ab' or p[-2:] == 'ba'):
                ints.append(0)
            else:
                ints.append(data_int[str(t)][p])
        data_dict[p] = ints

    to_fit = pd.DataFrame(data_dict, index = time)
    # print(to_fit)

    for i in auto_cross:
        to_fit[i] = to_fit[i].str.replace(',','').astype(float)

    data_dup = to_fit.copy(deep=True)

    for i in res:
        data_dup[i+'a'] = to_fit[i+'a']/to_fit[i+'a'][0.00]
        data_dup[i+'b'] = to_fit[i+'b']/to_fit[i+'b'][0.00]
        data_dup[i+'ba'] = to_fit[i+'ba']/to_fit[i+'a'][0.00]
        data_dup[i+'ab'] = to_fit[i+'ab']/to_fit[i+'b'][0.00]
    data_dup=data_dup.fillna(0)
    tag = file.split('_')[0]
    # print(data_dup)

    data_out = data_dup.copy(deep=True)
    # print(data_out)

    plt.style.use('ast-classic')
    # plt.rcParams.update({"font.family": "Times New Roman", "font.size":20})

    pdf = matplotlib.backends.backend_pdf.PdfPages('OsUBQ_zz_paper_'+tag+'.pdf')

    for i in range(len(auto_cross)):
        if i%4 == 0:
            comboY = np.array([])
            y = np.array(data_dup[auto_cross[i]], dtype='float')
            comboY = np.hstack((comboY, y))
        else:
            y = np.array(data_dup[auto_cross[i]], dtype='float')
            comboY = np.hstack((comboY, y))

        if (i != 0) and ((i+1)%4 == 0):
            
            initial_param = np.array([1,1,1,1])
            fittedparam, pcov = curve_fit(comboFunc, comboX, comboY, initial_param,method='lm')

            r1a,r1b,kab,kba = fittedparam
            perr = np.sqrt(np.diag(pcov))
            # value +/- std_err
            print("{r[0]:.2f}\u00B1{s[0]:0.2f}, {r[1]:0.2f}\u00B1{s[1]:0.2f}, {r[2]:0.2f}\u00B1{s[2]:0.2f}, {r[3]:0.2f}\u00B1{s[3]:0.2f}".format(r=fittedparam, s=perr))
            
            text = ('$R_{1a}$ = %.2f\u00B1%.2f $s^{-1}$ \n$R_{1b}$ = %.2f\u00B1%.2f $s^{-1}$ \n$k_{ab}$ = %.2f\u00B1%.2f $s^{-1}$ \n$k_{ba}$ = %.2f\u00B1%.2f $s^{-1}$'%(r1a, perr[0], r1b, perr[1], kab, perr[2], kba, perr[3]))

            # print(len(comboFunc(comboX, *fittedparam)))
            # print(len(comboY))
            # chi=np.sum(((comboFunc(comboX, *fittedparam) - comboY)/np.mean(perr))**2)
            # print('$CHI^{2}$ = %f'%(chi))

            fit_y1 = auto_aa(time_arr,r1a,r1b,kab,kba)
            fit_y2 = auto_bb(time_arr,r1a,r1b,kab,kba)
            fit_y3 = cross_ab(time_arr,r1a,r1b,kab,kba)
            fit_y4 = cross_ba(time_arr,r1a,r1b,kab,kba)

            # fig = plt.figure(figsize=(12,10), dpi=200)
            fig = plt.figure()

            units= "s\u207B\u00b9"
#            fig.suptitle('R1a = {r[0]:.2f}\u00B1{s[0]:0.2f} {g}, R1b = {r[1]:0.2f}\u00B1{s[1]:0.2f} {g}, kab = {r[2]:0.2f}\u00B1{s[2]:0.2f} {g}, kba = {r[3]:0.2f}\u00B1{s[3]:0.2f} {g}'.format(r=fittedparam, s=perr, g=units), y=0.95)

            res_id = auto_cross[i][:-2]

# write fit data to a DataFrame for outfile
            dd = np.vstack((fit_y1,fit_y2,fit_y4,fit_y3))
            col_name = ['F'+res_id+'a', 'F'+res_id+'b', 'F'+res_id+'ab', 'F'+res_id+'ba']
            fit_df = pd.DataFrame(dd.T, index=time, columns = col_name)
            # print(fit_df)

            data_out = pd.concat([data_out, fit_df], axis=1)

            fig.suptitle('%s(%s)'%(res_id, tag[:-1]+'oC'), x=0.5, y=0.93)

            plt.errorbar(comboX[:no],comboY[:no], yerr=abs(comboY[:no]-fit_y1), label='$%s_{a}$'%(res_id), c='black', fmt='o', markersize=10)
            # plt.scatter(comboX[:no],comboY[:no], label='$%s_{a}$'%(res_id), c='black', facecolors='none')
            plt.plot(time, fit_y1, 'black')
            plt.errorbar(comboX[no:2*no],comboY[no:2*no], yerr=abs(comboY[no:2*no]-fit_y2), label='$%s_{b}$'%(res_id), c='blue', fmt='o', markersize=10)
            # plt.scatter(comboX[no:2*no],comboY[no:2*no], label='$%s_{b}$'%(res_id), c='blue', facecolors='none')
            plt.plot(time, fit_y2, 'blue')
            plt.errorbar(comboX[2*no:3*no],comboY[2*no:3*no], yerr=abs(comboY[2*no:3*no]-fit_y3), label='$%s_{ab}$'%(res_id), c='red', fmt='s', markersize=10)
            # plt.scatter(comboX[2*no:3*no],comboY[2*no:3*no], label='$%s_{ab}$'%(res_id), c='red', marker='s', facecolors='none')
            plt.plot(time, fit_y3, 'red')
            plt.errorbar(comboX[3*no:4*no],comboY[3*no:4*no], yerr=abs(comboY[3*no:4*no]-fit_y4), label='$%s_{ba}$'%(res_id), c='green', fmt='s', markersize=10)
            # plt.scatter(comboX[3*no:4*no],comboY[3*no:4*no], label='$%s_{ba}$'%(res_id), c='green', marker='s', facecolors='none')
            plt.plot(time, fit_y4, 'green')
            plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
            plt.xlabel('$t (sec)$')
            plt.ylabel('$I_{t}$/$I_{0}$')
            plt.xlim(-0.05,1.05)
            plt.ylim(-0.05,1.05)
            # plt.legend(loc=1)
            plt.text(0.75,0.8,text, ha='left')
           
            plt.savefig(auto_cross[i]+'_paper_'+tag+'.svg', format='svg')
            pdf.savefig(fig)
            plt.clf()
            plt.close()

            # fo.write("{}\t{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}\n".format(auto_cross[i][:-2],*fittedparam))
            fo.write("{}\t{r[0]:.2f},{s[0]:0.2f}\t{r[1]:0.2f},{s[1]:0.2f}\t{r[2]:0.2f},{s[2]:0.2f}\t{r[3]:0.2f},{s[3]:0.2f}\n".format(auto_cross[i][:-2]+'_'+tag,r=fittedparam,s=perr))
    data_out = data_out.reset_index(level=0, names=['time'])
    # print(data_out.info())
    data_out.to_csv('master_'+tag+'.tab', sep='\t', index=False)
    pdf.close()
fo.close()
