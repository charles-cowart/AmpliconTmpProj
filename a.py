:52 (qiita@b2-001):AmpliconDevelopment$ cat p.py
import sys
import os
import re


miseq_root = '/sequencing/seqmount/KL_MiSeq_Runs'


def process_reads(reads):
        # extract all read elements as a list.
        # the contents of each Read element are highly regular.
        # for now, process w/out installing xml2dict or other
        # library into Qiita env.
        results = re.findall('<Read (.+?) \/>', reads)

        l = []
        for result in results:
            attributes = result.split(' ')
            d = {}
            for attribute in attributes:
                k, v = attribute.split('=')
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


def process_run_info_file(path):
    with open(path, 'r') as f:
        s = f.read()
        reads = re.search('<Reads>(.+?)</Reads>', s.replace('\n', ''))
        if reads:
            result = reads.group(1)
        else:
            raise ValueError("Cannot extract read information")
        return process_reads(result)


def main2(run_id):
    fp = os.path.join(miseq_root, run_id)
    if os.path.exists(fp):
        reads = process_run_info_file(os.path.join(fp, 'RunInfo.xml'))
        results = [x for x in reads if x['IsIndexedRead'] is True]

        if len(results) == 1:
            index1 = 'A' * results[0]['NumCycles']
            print("index1: %s" % index1)
        elif len(results) == 2:
            index1 = 'A' * results[0]['NumCycles']
            index2 = 'G' * results[1]['NumCycles']
            print("index1: %s" % index1)
            print("index2: %s" % index2)
        else:
            raise ValueError("too many indexed reads: %d" % len(results))
    else:
        raise ValueError("run_dir %s not found." % fp)




def main():
    for run_dir in os.listdir(miseq_root):
        try:
            print("processing %s" % run_dir)
            main2(run_dir)
            print("end processing %s" % run_dir)
            print("")
        except ValueError as e:
            print("Couldn't process %s: %s" % (run_dir, e))
main()
