import os, re
import pandas as pd


# Generic function to add suffix to the 'cell' column in DataFrame.
# suffix is equivalent to codename
def add_suffix(df, suffix):
    df.cell = df.cell + suffix


if __name__ == '__main__':

    # Set up some default parameters, i.e. working directory and filename
    parent_wd = '/media/luolab/ZA1BT1ER/yanting/vM21/'
    data_wd = '/media/luolab/ZA1BT1ER/yanting/vM21/mapping/'

    # Change working directory
    os.chdir(parent_wd)

    # Import experimental design, which stores library name, index, codename & experiment setup.
    exp_design = pd.read_excel('experimental_design_all.xlsx')

    # Iteratively enter each library
    for folder in os.listdir(data_wd):

        library_wd = os.path.join(data_wd, folder)
        os.chdir(library_wd)

        print('parsing ' + folder + ' ...')

        # Extract names from folder
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', folder)

        # Import data
        whitelist = pd.read_csv('_'.join([match.group(1), 'whitelist_washed.txt']), sep="\t",
                                names=['cell', 'candidate', 'Nreads', 'Ncandidate'])
        whitelist_Nuniquemap = pd.read_csv('_'.join([match.group(1), 'Nuniqmapped.txt']),
                                           delim_whitespace=True, names=['Nuniquemap', 'cell'])
        # CaiT pipeline uses 'id' as header
        # whitelist = pd.read_csv('whitelist.txt', sep="\t", names=['id', 'candidate', 'Nreads', 'Ncandidate'])
        # whitelist_Nuniquemap = pd.read_csv('whitelist.txt.mapped', delim_whitespace=True, names=['Nuniquemap', 'id'])

        # Find codename w.r.t. lib_prefix to use as suffix
        codename = exp_design.codename[exp_design.lib_prefix == match.group(1)]

        # Use iloc to access by position rather than label (awesome!)
        suffix = '_' + codename.iloc[0]

        # Add suffix
        add_suffix(whitelist, suffix)
        add_suffix(whitelist_Nuniquemap, suffix)

        # Construct mapping stat DataFrame
        frames = [whitelist.Nreads, whitelist_Nuniquemap.Nuniquemap, whitelist.cell]
        mapping_stat = pd.concat(frames, axis=1)

        # Output
        out_nameparts = [match.group(1), 'mapping_stats.txt']
        out_filename = '_'.join(out_nameparts)
        mapping_stat.to_csv(out_filename, sep="\t", index=False)
