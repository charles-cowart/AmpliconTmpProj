import sample_sheet
import pandas as pd
from metapool import KLSampleSheet


def generate_amplicon_sample_sheet():
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
    sheet.Reads = [151, 151]
    sheet.Settings['OverrideCycles'] = 'Y151;I12;Y151'

    sheet.Settings['MaskShortReads'] = '1'

    dummy_samples = {'Sample_ID': 'Bacterial_Test_Study_Plates1_2',
                     'Sample_Name': '',
                     'Sample_Plate': '',
                     'Sample_Well': '',
                     'I7_Index_ID': '',
                     'index': 'AAAAAAAA',
                     'I5_Index_ID': '',
                     'index2': 'GGGGGGGG',
                     'Sample_Project': 'Bacterial_Test_Study_Plates1_2',
                     'Description': ''
                     }
    sheet.add_sample(sample_sheet.Sample(dummy_samples))

    # we'll get these from input parameters as well.
    contacts = [['test@lol.com', 'Baz'], ['tester@rofl.com', 'FooBar_666']]
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

    with open('out.csv', 'w') as f:
        sheet.write(f, 1)

