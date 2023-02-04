import sys
import os
import re
import sample_sheet
import pandas as pd
from metapool import KLSampleSheet


miseq_root = '/sequencing/seqmount/KL_MiSeq_Runs'


def generate_amplicon_sample_sheet(read1, read2, override_cycles, index1,
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
    sheet.Reads = [read1, read2]
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

        print(l)
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
        reads = process_run_info_file(os.path.join(fp, 'RunInfo.xml'))
        results = [x for x in reads if x['IsIndexedRead'] is True]

        if len(results) == 1:
            index1 = 'A' * results[0]['NumCycles']
            index2 = None
        elif len(results) == 2:
            index1 = 'A' * results[0]['NumCycles']
            index2 = 'G' * results[1]['NumCycles']
        else:
            raise ValueError("too many indexed reads: %d" % len(results))
    else:
        raise ValueError("run_dir %s not found." % fp)

    return index1, index2


def process_run_dir(run_id, output_dir):
    index1, index2 = get_runinfo_params(run_id)

    # these will need to be pulled from runinfo.xml as well.
    read1 = 151
    read2 = 151
    override_cycles = 'Y151;I12;Y151'
    contacts = [['test@lol.com', 'Baz'], ['tester@rofl.com', 'FooBar_666']]
    sheet = generate_amplicon_sample_sheet(read1, read2, override_cycles,
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
            print("Couldn't process %s: %s" % (run_dir, e))
            results = 1


    return results


if __name__ == '__main__':
    # main will return a proper return code
    # sys.exit(main())
    process_run_dir('180404_A00169_0085_AH7CV2DMXX', '.')
