import matplotlib
import pandas as pd
try:
    matplotlib.use("Agg")
except ValueError:
    pass

import numpy, sys, pylab, os.path, pandas, seaborn
from tqdm import tqdm

WINDOW_SIZE = 20000

def parse_bed(bed,window=WINDOW_SIZE):
    chromos = bed['chr'].unique()
    for chromo in chromos:
        df = bed[bed['chr']==chromo] 
        df = df[['base','depth']]
        df['rolling']=df['depth'].rolling(window=window).mean()
        df['rolling']=df['rolling'].fillna(0)
        length = df.shape[0]
        yield chromo, length, df
        
def format_chr(chromo_id):
    if chromo_id.startswith('Ld'):
        chromo_id = chromo_id.split('_')[0]
        chromo_id = chromo_id.replace('Ld','')
        chromo_id = chromo_id.replace('kinetoplast','Kin')
        return chromo_id

def make_plot(bed, save_to):
    seaborn.set_style("white")
    df = pd.read_csv(bed, sep='\t',header=None,index_col=None)
    df.columns = ['chr', 'base', 'depth']
    offset = 0
    medians = {}
    with seaborn.color_palette("husl", len(df['chr'].unique())):
        for chromo, length, counts in tqdm(parse_bed(df), total=len(df['chr'].unique())):
            counts=counts.fillna(0)
            if "Un" in chromo or "random" in chromo or "hap" in chromo or "M" in chromo or "mt" in chromo:
                continue # skip weird chromosomes
            pylab.plot(counts['base']+offset, counts['rolling'])
            y_text = counts['rolling'].max()
            pylab.text(offset + length*0.5, y_text, format_chr(chromo), 
                       horizontalalignment="center", verticalalignment="center")
            offset += length

            medians[chromo] = counts['rolling'].median()
            pylab.axvline(x=offset, ymin=0, ymax=1, ls='-.', alpha = 0.5)

    #pylab.legend()
    pylab.xlabel("bp")
    pylab.ylabel("coverage (%d-bp averaged)" % WINDOW_SIZE)

    pylab.gca().set_xlim([0-WINDOW_SIZE*10, offset+WINDOW_SIZE*10])
    pylab.title(bed)
    pylab.gcf().set_size_inches(16, 4)
    pylab.savefig(save_to)

if __name__ == '__main__':
    bed = sys.argv[1]
    save_to = sys.argv[2]
    make_plot(bed,save_to)
