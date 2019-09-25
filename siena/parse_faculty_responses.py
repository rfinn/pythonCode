#!/usr/bin/env python

import csv
import argparse
from pylab import * 
import os

def make_hist(x,figname,pltitle,plotsingle=False):
    if plotsingle:
        figure()
    mybins=arange(0,2.7,.5)-.25
    #print mybins
    #print x
    t=hist(x,bins=mybins)
    xlim(-0.5,2.5)
    ylim(0,30)
    xticks(arange(3),('$Agree$','$Disagree$','$Undecided$'),fontsize=12)
    title(pltitle,fontsize=14)
    text(.95,.9,'$N_{vote} = %i$'%(sum(t[0])),transform=gca().transAxes,horizontalalignment='right')
    if plotsingle:
        savefig(figname)
def make_time_hist(x,figname,pltitle,xtick_names,plotsingle=False):
    if plotsingle:
        figure(figsize=(8,4))
    mybins=arange(len(x))
    #print mybins
    #print x
    bar(mybins,x)
    print t
    #xlim(-0.5,2.5)
    ylim(0,17)

    xtick_names=['Mon PM', ' Tues PM', ' Wed PM', ' Thurs PM', ' Fri PM', ' Sat PM', ' Sun PM',' Sat',' Sun' ]
    xticks(arange(9)+.4,(xtick_names),fontsize=8)
    yticks(arange(0,17,5))
    title(pltitle,fontsize=14)
    text(.95,.9,'$N_{vote} = %i$'%(sum(x)),transform=gca().transAxes,horizontalalignment='right')
    if plotsingle:
        savefig(figname)

class texfile():
    def __init__(self,filename):
        self.filename=filename
        self.outfile=open(filename,'w')
        self.outfile.write('\documentclass[11pt]{article} \n')
        self.outfile.write('\usepackage{graphicx} \n')
        self.outfile.write('\setlength{\\textwidth}{6.5in} \n')
        self.outfile.write('\setlength{\\textheight}{9in} \n')
        self.outfile.write('\setlength{\\topmargin}{-1.cm} \n')
        self.outfile.write('\setlength{\\oddsidemargin}{0in} \n')
        self.outfile.write('\setlength{\\evensidemargin}{0in} \n')        
        self.outfile.write('\usepackage{graphicx} \n')
        self.outfile.write('\\begin{document} \n')

    def add_section(self,title):
        self.outfile.write('\section{'+title+'} \n')
    def start_list(self):
        self.outfile.write('\\begin{itemize} \n')
    def add_item(self,item):
        self.outfile.write('\item '+item+'\n')
    def end_list(self):
        self.outfile.write('\\end{itemize} \n')
    def add_figure(self,figname,caption=None):
        self.outfile.write('\\begin{figure}[h] \n')
        self.outfile.write('\includegraphics[width=\\textwidth]{'+figname+'} \n')
        if caption:
            self.outfile.write('\caption{'+caption+'} \n')
        self.outfile.write('\end{figure} \n')

    def close_file(self):
        self.outfile.write('\end{document} \n')
        self.outfile.close()
    def makepdf(self):
        ftrunk=self.filename.split('.')[0]
        os.system('xelatex '+ftrunk)

################################################################################
if __name__ == "__main__":

    ############################################################################
    # Parse the arguments
    ############################################################################
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('input_file_name', type=str, default=None, 
            help='Input file name')
    '''
    parser.add_argument('--nofinal', dest='nofinal', action='store_true',
            default=True,
            help='Dump the names of the students and an index for each.')
    parser.add_argument('--dump-names', dest='dump_names', action='store_true',
            default=True,
            help='Dump the names of the students and an index for each.')
    '''
    args = parser.parse_args()
    filename = args.input_file_name
    t=open(filename, 'rb')
    infile = csv.reader(t,delimiter='{', quotechar='#')

    senate = []
    alt_time = []
    best_times=zeros(9,'f')
    '''
    0 = Mon Eve
    1 = Tues Eve
    2 = Wed Eve
    3 = Thurs Eve
    4 = Fri Eve
    5 = Sat Eve
    6 = Sun even
    7 = sat
    8 = Sun 
    '''
    best_time_names=['Monday evening', ' Tuesday evening', ' Wednesday evening', ' Thursday evening', ' Friday evening', ' Saturday evening', ' Sunday evening',' Saturday',' Sunday' ]

    line_num = 0
    comments=texfile('comments.tex')
    comments.add_section('Survey Results')
    comments.add_figure('responses.png',caption='Survey Results for (left) finding an alternate meeting time and (right) forming a faculty senate.')
    comments.add_figure('best_times.png',caption='Preferred times for holding a meeting during the evening or weekend.')

    comments.add_section('Other Suggestions')
    comments.start_list()
    for row in infile:

        if line_num==0:
            line_num += 1
            continue
        #print len(row),row
        try:
            if row[1].find('Agree') > -1:
                alt_time.append(0)
            elif row[1].find('Disagree') > -1:
                alt_time.append(1)
            elif row[1].find('Undecided') > -1:
                alt_time.append(2)
            if row[2].find('Agree') > -1:
                senate.append(0)
            elif row[2].find('Disagree') > -1:
                senate.append(1)
            elif row[2].find('Undecided') > -1:
                senate.append(2)
            # add comments to a tex file
            if len(row[3]) > 2:
                print 'writing comment'
                comments.add_item(str(row[3]))
            if len(row[4]) > 1:
                t=row[4].split(',')
                #print t
                for d in t:
                    i=0
                    for n in best_time_names:
                        print d, n
                        if d == n:#(d.find(n) > -1):
                            best_times[i] += 1
                            print 'found match with ',n
                            break
                        i += 1
                            
        except IndexError:
            continue
        line_num += 1
alt_time=array(alt_time,'i')
senate=array(senate,'i')
figure(figsize=(8,3))
subplots_adjust(wspace=.4,left=.1,right=.9)
subplot(1,2,1)
make_hist(alt_time,'alt_time.png','$Alternate \ Meeting \ Time$',plotsingle=False)
subplot(1,2,2)
make_hist(senate,'senate.png','$Form \ a \ Faculty \ Senate$',plotsingle=False)
savefig('responses.png')

make_time_hist(best_times,'best_times.png','$Best \ Times \ for \ Meeting$', best_time_names,plotsingle=True)
comments.end_list()
comments.close_file()
comments.makepdf()
