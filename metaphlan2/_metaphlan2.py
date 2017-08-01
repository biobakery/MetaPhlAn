import subprocess as sb
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.per_sample_sequences import SingleLanePerSampleSingleEndFastqDirFmt, SingleLanePerSamplePairedEndFastqDirFmt
import tempfile
import biom

def metaphlan2_helper(raw_data, nproc, input_type, output_file, verbose=True):
    cmd = ['metaphlan2.py', '--input_type', str(input_type), '--biom',
           '--nproc', str(nproc), str(raw_data)]

    if verbose:
        print("Running external command line application. This may print messages to stdout and/or stderr.", end='\n\n')
        print("Command: {} > {}".format(' '.join(cmd), output_file), end='\n\n')

    with open(str(output_file), 'w') as output_file_h:
        sb.run(cmd, stdout=output_file_h, check=True)


def profile_single_fastq(raw_data: SingleLanePerSampleSingleEndFastqDirFmt, nproc: int=1) -> biom.Table:
    with tempfile.TemporaryDirectory() as temp_dir_name:
        biom_fp = os.path.join(temp_dir_name, 'output.tsv.biom')

        metaphlan2_helper(raw_data, nproc, input_type='multifastq', output_file=biom_fp)
        table = biom.load_table(biom_fp)
    return table


def profile_paired_fastq(raw_data: SingleLanePerSamplePairedEndFastqDirFmt, nproc: int=1) -> biom.Table:
    with tempfile.TemporaryDirectory() as temp_dir_name:
        biom_fp = os.path.join(temp_dir_name, 'output.tsv.biom')

        metaphlan2_helper(raw_data, nproc, input_type='multifastq', output_file=biom_fp)
        table = biom.load_table(biom_fp)
    return table


#def metaphlan2_single_fasta(raw_data: SingleLanePerSampleSingleEndFastaDirFmt, nproc: int=1) -> FeatureTable[Frequency]:
#    result = FeatureTable[Frequency]()
#    metaphlan2_helper(raw_data, nproc, 'multifasta', result)
#    return result


#def metaphlan2_paired_fasta(raw_data: SingleLanePerSamplePairedEndFastaDirFmt, nproc: int=1) -> FeatureTable[Frequency]:
#    result = FeatureTable[Frequency]()
#    metaphlan2_helper(raw_data, nproc, 'multifasta', result)
#    return result

