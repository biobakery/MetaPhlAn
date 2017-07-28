import subprocess as sb
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.per_sample_sequences import SingleLanePerSampleSingleEndFastqDirFmt, SingleLanePerSamplePairedEndFastqDirFmt


def metaphlan2_helper(raw_data, nproc, input_type, output_file, verbose=True):
    cmd = ['metaphlan2.py', '--input_type', str(input_type), '--nproc', str(nproc), str(raw_data)]

    if verbose:
        print("Running external command line application. This may print messages to stdout and/or stderr.", end='\n\n')
        print("Command: {} > {}".format(' '.join(cmd), output_file), end='\n\n')

    with open(str(output_file), 'w') as output_file_h:
        sb.run(cmd, stdout=output_file_h, check=True)


def metaphlan2_single_fastq(raw_data: SingleLanePerSampleSingleEndFastqDirFmt, nproc: int=1) -> FeatureTable[Frequency]:
    result = FeatureTable[Frequency]()
    metaphlan2_helper(raw_data, nproc, 'multifastq', result)
    return result


def metaphlan2_paired_fastq(raw_data: SingleLanePerSamplePairedEndFastqDirFmt, nproc: int=1) -> FeatureTable[Frequency]:
    result = FeatureTable[Frequency]()
    metaphlan2_helper(raw_data, nproc, 'multifastq', result)
    return result


#def metaphlan2_single_fasta(raw_data: SingleLanePerSampleSingleEndFastaDirFmt, nproc: int=1) -> FeatureTable[Frequency]:
#    result = FeatureTable[Frequency]()
#    metaphlan2_helper(raw_data, nproc, 'multifasta', result)
#    return result


#def metaphlan2_paired_fasta(raw_data: SingleLanePerSamplePairedEndFastaDirFmt, nproc: int=1) -> FeatureTable[Frequency]:
#    result = FeatureTable[Frequency]()
#    metaphlan2_helper(raw_data, nproc, 'multifasta', result)
#    return result

