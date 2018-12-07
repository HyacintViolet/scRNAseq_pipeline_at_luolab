import os
import re
import shutil


data_wd = '/media/luolab/ZA1BT1ER/yanting/vM4/vM4_CaiT/'
folders_list = os.listdir(data_wd)


# ----------------------------------------------------------------------------------------------------------------------
# Rename file name prefix: yantingx_yyyyyy.zzz -> YTxxxxxx_yyyyyyy.zzz
# ----------------------------------------------------------------------------------------------------------------------


for folder in folders_list:

    print('Parsing ' + folder)

    os.chdir(os.path.join(data_wd, folder))

    # Grab folder name
    match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', folder)

    files = os.listdir(os.path.join(data_wd, folder))

    for file in files:

        if file.startswith('yanting'):

            suffix = re.search('^(?:yanting)(?:.*)(?:[\d]+)_(.*)$', file)

            os.rename(file, '_'.join([match.group(1), suffix.group(1)]))

# ----------------------------------------------------------------------------------------------------------------------
# Recursively remove "_suffix_renamed.tsv", or "whitelistx.txt" and "Nuniqmapped.txt"
# ----------------------------------------------------------------------------------------------------------------------


# def remove_renamed(src):
#     for item in os.listdir(src):
#         s = os.path.join(src, item)
#         if os.path.isdir(s):
#             remove_renamed(s)
#         elif os.path.basename(s).endswith('suffix_renamed.tsv'):
#             os.remove(s)
#
#
# remove_renamed(data_wd)


# def remove_whitelist(src):
#     for item in os.listdir(src):
#         s = os.path.join(src, item)
#         if os.path.isdir(s):
#             remove_whitelist(s)
#         elif (os.path.basename(s) == 'whitelist.txt') or (os.path.basename(s) == 'Nuniqmapped.txt'):
#             os.remove(s)
#
#
# remove_whitelist(data_wd)

# ----------------------------------------------------------------------------------------------------------------------
# Duplicate whitelistx.txt -> whitelist.txt
# ----------------------------------------------------------------------------------------------------------------------


def navigate_and_rename_whitelist(src):
    for item in os.listdir(src):
        s = os.path.join(src, item)
        if os.path.isdir(s):
            navigate_and_rename_whitelist(s)
        elif os.path.basename(s).startswith('whitelist') and (not os.path.basename(s).endswith('mapped')):
            match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', os.path.basename(os.path.dirname(s)))
            print('Parsing' + s)
            shutil.copy(s, os.path.join(src, '_'.join([match.group(1), 'whitelist.txt'])))


navigate_and_rename_whitelist(data_wd)

# ----------------------------------------------------------------------------------------------------------------------
# Duplicate whitelistx.txt.mapped -> Nuniqmapped.txt
# ----------------------------------------------------------------------------------------------------------------------


def navigate_and_rename_nuniqmapped(src):
    for item in os.listdir(src):
        s = os.path.join(src, item)
        if os.path.isdir(s):
            navigate_and_rename_nuniqmapped(s)
        elif os.path.basename(s).endswith('mapped'):
            match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', os.path.basename(os.path.dirname(s)))
            print('Parsing' + s)
            shutil.copy(s, os.path.join(src, '_'.join([match.group(1), 'Nuniqmapped.txt'])))


navigate_and_rename_nuniqmapped(data_wd)
