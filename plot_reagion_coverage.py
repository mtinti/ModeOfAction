import sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import gc
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
plt.style.use('ggplot')
from adjustText import adjust_text
import numpy as np

def make_dict(cell):
    temp = {}
    for item in cell.split(';'):
        temp[item.split('=')[0].strip()]=item.split('=')[1].strip()
    return temp

def get_field(cell, key):
    cell_dict = make_dict(cell)
    return cell_dict[key].replace('%2','')
        
def gff_to_pandas(infile):
    df = pd.read_csv(infile, sep='\t', header=None, index_col=None,comment='#')
    df.columns = ['gene_chr','source','ftype','gene_start',
                  'gene_end', 'score', 'gene_strand', 'score2', 'gene_desc']
    df=df[df['ftype']=='gene']
    df['gene_id']=[get_field(cell,key='ID') for cell in df['gene_desc']]
    df['gene_desc']=[get_field(cell,key='description') for cell in df['gene_desc']]
    df.drop(['score','score2','source','ftype'],axis=1,inplace=True)
    return df


def add_track(bed_file = str, chro = str, 
              start =int, end=int, ax=None,
              label = str, color=str, roll=False):
    df = pd.read_csv(bed_file, sep='\t', header=None,index_col=None)
    if bed_file.endswith('_d.bed'):
        df.columns = ['chr', 'base', 'coverage']
    if bed_file.endswith('_bg.bed'):
        df.columns = ['chr', 'base', 'to_base', 'coverage']
        
    test = df[df['chr']==chro]
    del df
    gc.collect()
    if roll:
        test['coverage']=test['coverage']-test['coverage'].mean()
        test['coverage']=test['coverage'].rolling(window=2000).median()
        #test=test[test['coverage']>0]

    test=test[ (test['base']>start) & (test['base']<end)]    
    test.plot(x='base', y='coverage', ax=ax, label=label, c=color)
    return test.base.min(), test.base.max()

def plot_region(coverage_file, ff_file, fr_file, rf_file, rr_file, chro, start, end, save_to, counts, peak_index, gff):

    f = plt.figure(figsize=(16,6))
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    ax = plt.subplot(gs[0])
    axt = ax.twinx()

    _, _, = add_track(bed_file = coverage_file, chro = chro, 
              start =start, end=end, ax=ax, label='coverage',color='b', roll=True)


    ax.set_ylabel('Coverage')
    axt.set_ylabel('Primer Coverage')


    region = gff[(gff['gene_chr']==chro) & (gff['gene_start'] > start ) & (gff['gene_end'] < end )]

    ax.set_xlim(start, end)
    ax.set_ylim(bottom=0)
    
    ax2 = plt.subplot(gs[1], sharex = ax)

    add_track(bed_file = ff_file, chro = chro, 
              start =start, end=end, ax=axt, label='ff',color='r')

    add_track(bed_file = fr_file, chro = chro, 
              start =start, end=end, ax=axt, label='fr',color='g')


    add_track(bed_file = rf_file, chro = chro, 
              start =start, end=end, ax=axt, label='rf',color='y')

    add_track(bed_file = rr_file, chro = chro, 
              start =start, end=end, ax=axt, label='rr',color='k')


    ax.legend(loc="upper left")


    #print(region.head())
    for gene in region.index.values:
        #print(gene,'__\n\n')
        gid = region.loc[gene].name
        gstart = region.loc[gene].gene_start
        gend = region.loc[gene].gene_end
        ax2.add_patch(Rectangle((gstart, 0), width=gend-gstart, height=0.1))

    texts = []
    for i in region.index.values:
        #print(i)
        if i in counts.index:
            texts.append(ax2.text(region.loc[i].gene_start, 0.1, i))
                            #if i in counts.index ]  

    x= np.arange(start, end, 100) 
    y=[0.2 for n in x]
    ax2.scatter(x, y,alpha=0)
    ax2.set_ylim(0,1)
    adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'), ax=ax2)
    title = 'peak {peak_index} loc: {chro}:{start}-{end}'.format(chro=chro,start=start,end=end,peak_index=peak_index)
    ax.set_title(title)
    plt.savefig(save_to+title+'.png')
    

if __name__ == '__main__':

    gff = sys.argv[1]
    gff = gff_to_pandas(gff)
    gff.head()
    
    peaks = pd.read_table(sys.argv[2], sep='\t', comment='#')
    peaks = peaks.sort_values('pileup', ascending=False)

    counts = pd.read_table(sys.argv[9], sep='\t', index_col=[0])
    counts = counts[counts[counts.columns[0]]>0]
    counts = counts[counts[counts.columns[0]] >counts[counts.columns[0]].mean() + counts[counts.columns[0]].std()]
    
    coverage_file = sys.argv[3]
    ff_file = sys.argv[4]
    fr_file = sys.argv[5]
    rf_file = sys.argv[6]
    rr_file = sys.argv[7]
    save_to = sys.argv[8]
    for peak_index, peak in enumerate(peaks.index.values[0:3]):
        chro = peaks.loc[peak]['chr']
        start = peaks.loc[peak]['start']
        end = peaks.loc[peak]['end']
        plot_region(coverage_file, ff_file, fr_file, rf_file, rr_file, chro, start, end, save_to, counts, peak_index+1, gff)
