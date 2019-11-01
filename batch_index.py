cmd_index = []
    for out in os.listdir(dst):

        # Setup output directory
        out_dir = os.path.join(dst, out)

        # Change directory to output folder
        os.chdir(out_dir)

        # Grab folder name
        match = re.search('^([^_]*)_([^_]*)_([^_]*)_([^_]*)$', out)out_name_samtools
        prefix = match.group(1)

        # Input file name for samtools index
        out_name_samtools = '_'.join([prefix, 'assigned_sorted.bam'])
        index_in = os.path.join(out_dir, out_name_samtools)

        # Construct command
        cmd_this_index = 'samtools index ' + index_in
        cmd_index.append(cmd_this_index)

    # Parallel run by Pool
    pool = mp.Pool(16)
    pool.map(work, cmd_index)
    print('samtools index: finished.')