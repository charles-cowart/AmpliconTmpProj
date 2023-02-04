import sys
import os
import re
import sample_sheet
import pandas as pd
from metapool import KLSampleSheet


miseq_root = '/sequencing/seqmount/KL_MiSeq_Runs'


def generate_amplicon_sample_sheet(reads, override_cycles, index1,
                                   index2, contacts):
    # create object and initialize header
    sheet = KLSampleSheet()
    sheet.Header['IEMFileVersion'] = '4'
    sheet.Header['Date'] = '10/27/22'
    sheet.Header['Workflow'] = 'GenerateFASTQ'
    sheet.Header['Application'] = 'FASTQ Only'
    sheet.Header['Assay'] = 'TruSeq HT'
    sheet.Header['Description'] = 'test_run'
    sheet.Header['Chemistry'] = 'Amplicon'

    # set Reads and Settings according to input values
    # we'll get this from the code on the server
    sheet.Reads = reads
    sheet.Settings['OverrideCycles'] = override_cycles

    sheet.Settings['MaskShortReads'] = '1'

    if index1 is None:
        index1 = ''

    if index2 is None:
        index2 = ''

    dummy_samples = {'Sample_ID': 'Bacterial_Test_Study_Plates1_2',
                     'Sample_Name': '',
                     'Sample_Plate': '',
                     'Sample_Well': '',
                     'I7_Index_ID': '',
                     'index': index1,
                     'I5_Index_ID': '',
                     'index2': index2,
                     'Sample_Project': 'Bacterial_Test_Study_Plates1_2',
                     'Description': ''
                     }
    sheet.add_sample(sample_sheet.Sample(dummy_samples))

    # we'll get these from input parameters as well.
    contacts = pd.DataFrame(columns=['Email', 'Sample_Project'],
                            data=contacts)
    sheet.Contact = contacts

    # add a dummy sample.
    samples = [['Bacterial_Test_Study_Plates1_2', 'NA', 'NA',
                'FALSE', 'FALSE', '14782']]

    samples = pd.DataFrame(columns=['Project', 'ForwardAdapter',
                                    'ReverseAdapter', 'PolyGTrimming',
                                    'HumanFiltering', 'QiitaID'],
                           data=samples)

    sheet.Bioinformatics = samples

    return sheet


def process_run_info_file(run_info_fp):
    def process_reads(reads):
        # extract all read elements as a list.
        # the contents of each Read element are highly regular.
        # for now, process w/out installing xml2dict or other
        # library into Qiita env.
        reads = reads.split('<Read ')
        reads = [x.strip() for x in reads]
        reads = [x for x in reads if x != '']

        l = []
        for read in reads:
            d = {}
            read = read.strip('<').strip('>').strip('/')
            pairs = read.split(' ')
            pairs = [x for x in pairs if x != '']
            for pair in pairs:
                k, v = pair.split('=')
                if k in ['NumCycles', 'Number']:
                    v = int(v.strip('"'))
                elif k in ['IsIndexedRead']:
                    v = v.strip('"')
                    v = False if v == 'N' else True
                else:
                    raise ValueError("Unknown key: %s" % k)

                d[k] = v
            l.append(d)

        return l

    with open(run_info_fp, 'r') as f:
        s = f.read()
        s = s.replace('\n', '')
        reads = re.search('<Reads>(.+?)</Reads>', s)
        if reads:
            result = reads.group(1)
        else:
            raise ValueError("Cannot extract read information")
        result = process_reads(result)
        return result


def get_runinfo_params(run_id):
    fp = os.path.join(miseq_root, run_id)
    if os.path.exists(fp):
        print(os.path.join(fp, 'RunInfo.xml'))
        reads = process_run_info_file(os.path.join(fp, 'RunInfo.xml'))

        reads_list = [x['NumCycles'] for x in reads if x['IsIndexedRead'] is False]

        if len(reads_list) != 2:
            raise ValueError("RunInfo.xml contains an unexpected number of Non-indexed Reads (%d)" % len(reads_list))

        results = [x for x in reads if x['IsIndexedRead'] is True]

        if len(results) == 1:
            index1 = 'A' * results[0]['NumCycles']
            index2 = None
        elif len(results) == 2:
            index1 = 'A' * results[0]['NumCycles']
            index2 = 'G' * results[1]['NumCycles']
        elif len(results) == 0:
            raise ValueError("RunInfo.xml does not contain indexed reads")
        else:
            raise ValueError("too many indexed reads: %d" % len(results))
    else:
        raise ValueError("run_dir %s not found." % fp)

    return index1, index2, reads_list


def process_run_dir(run_id, output_dir):
    index1, index2, reads = get_runinfo_params(run_id)

    if index2:
        override_cycles = f'Y{reads[0]};I{len(index1)};I{len(index2)};Y{reads[1]}'
    else:
        override_cycles = f'Y{reads[0]};I{len(index1)};Y{reads[1]}'

    contacts = [['test@lol.com', 'Baz'], ['tester@rofl.com', 'FooBar_666']]
    sheet = generate_amplicon_sample_sheet(reads, override_cycles,
                                           index1, index2, contacts)

    result_path = os.path.join(output_dir, run_id + '.csv')

    with open(result_path, 'w') as f:
        sheet.write(f, 1)

    return result_path


def main():
    results = 0
    for run_dir in os.listdir(miseq_root):
        if run_dir in ['.DS_Store', '170425_M05314_0006_000000000-AEARP', '170915_M05314_0028_000000000-BDMKR']:
            # 170ddd directories don't have a runinfo.xml file.
            continue

        try:
            print("processing %s" % run_dir)
            sample_sheet_fp = process_run_dir(run_dir, './test_output_files')
            print("Created %s" % sample_sheet_fp)
        except (FileNotFoundError, ValueError) as e:
            # can also be because of a single non-indexed read when two are expected.
            print("Couldn't process %s: %s (or permissions not set)" % (run_dir, e))
            results = 1


    return results


if __name__ == '__main__':
    # main will return a proper return code
    sys.exit(main())


