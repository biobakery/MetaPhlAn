import subprocess as sb
from q2_types.per_sample_sequences import SingleLanePerSampleSingleEndFastqDirFmt
from q2_types.per_sample_sequences import SingleLanePerSamplePairedEndFastqDirFmt
import tempfile
import biom
import os
# import glob
# import gzip


def metaphlan2_helper(raw_data, nproc, input_type, output_file, verbose=True):
    # cmd = ['metaphlan2.py', '--input_type', str(input_type), '--biom',
    #        str(output_file), '--nproc', str(nproc)]
    cmd = ['metaphlan2.py', str(raw_data), '--input_type', str(input_type),
           '--biom', str(output_file), '--nproc', str(nproc)]

    if verbose:
        print("\nRunning external command line application. This may print "
              "messages to stdout and/or stderr.")
        print("Command: {}".format(' '.join(cmd)), end='\n\n')

    # for f in glob.iglob(os.path.join(str(raw_data), '*.fastq.gz')):
    #     print('Processing "{}"...'.format(f))
    #     t = sb.Popen(cmd, stdin=sb.PIPE, stdout=sb.DEVNULL)
    #     t.communicate(gzip.open(f, "rb").read())
    sb.run(cmd, check=True)


def profile_single_fastq(raw_data: SingleLanePerSampleSingleEndFastqDirFmt,
                         nproc: int=1) -> biom.Table:
    output_biom = None

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_output_biom = os.path.join(tmp_dir, 'mp2_tmp_output.biom')
        metaphlan2_helper(raw_data, nproc, 'multifastq', tmp_output_biom)
        output_biom = biom.load_table(tmp_output_biom)

    return output_biom


def profile_paired_fastq(raw_data: SingleLanePerSamplePairedEndFastqDirFmt,
                         nproc: int=1) -> biom.Table:
    output_biom = None

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_output_biom = os.path.join(tmp_dir, 'mp2_tmp_output.biom')
        metaphlan2_helper(raw_data, nproc, 'multifastq', tmp_output_biom)
        output_biom = biom.load_table(tmp_output_biom)

    return output_biom
