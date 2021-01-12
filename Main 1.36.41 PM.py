import os

import pandas as pd
from collections import OrderedDict
import gzip
import numpy as np
import sys

# print(list(data.columns.values))

VCF_HEADER = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

comments = 0

filename = sys.argv[1]
output_file = sys.argv[2]

while filename and output_file is None | filename or output_file is not None:
    print("Please provide valid path to the input file as a first argument and output file as a second argument")
    filename = sys.argv[1]
    output_file = sys.argv[2]

if os.path.exists(filename) is not True:
    print("The provided path to the input file is not valid. Please provide valid path to the input file.")

def dataframe(filename, large=True):
    """Open an optionally gzipped VCF file and return a pandas.DataFrame with
    each INFO field included as a column in the dataframe.
    Note: Using large=False with large VCF files. It will be painfully slow.
    :param filename:    An optionally gzipped VCF file.
    :param large:       Use this with large VCF files to skip the ## lines and
                        leave the INFO fields unseparated as a single column.
    """
    if large:
        # Set the proper argument if the file is compressed.
        comp = 'gzip' \
            if filename.endswith('.gz') \
            else None
        # Count how many comment lines should be skipped.
        comments = _count_comments(filename)
        # Return a simple DataFrame without splitting the INFO column.
        return pd.read_table(filename, compression=comp, skiprows=comments, names=VCF_HEADER, usecols=range(7))

    # Each column is a list stored as a value in this dict. The keys for this
    # dict are the VCF column names and the keys in the INFO column.
    result = OrderedDict()
    # Parse each line in the VCF file into a dict.
    for i, line in enumerate(lines(filename)):
        # with open(filename) as fh:
        # for i, line in fh:
        for key in line.keys():
            # This key has not been seen yet, so set it to None for all
            # previous lines.
            if key not in result:
                result[key] = [None] * i
        # Ensure this row has some value for each column.
        for key in result.keys():
            result[key].append(line.get(key, None))

    test = pd.DataFrame(result)
    return test


def add_samples(filename, info_data):
    comments = _count_comments(filename) - 1
    "This function will add the samples to the parsed dataframe"
    data_1 = len(pd.read_csv(filename, sep='\t', header=comments))
    data_2 = data_1[data_1.columns[8:]]

    final_data = info_data.join(data_2)
    return final_data


def lines(filename):
    """Open an optionally gzipped VCF file and generate an OrderedDict for
    each line.
    """
    fn_open = gzip.open if filename.endswith('.gz') else open

    with fn_open(filename) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            else:
                yield parse(line)


def parse(line):
    """Parse a single VCF line and return an OrderedDict.
    """
    result = OrderedDict()

    fields = line.rstrip().split('\t')

    # Read the values in the first seven columns.
    for i, col in enumerate(VCF_HEADER[:7]):
        result[col] = _get_value(fields[i])

    # INFO field consists of "key1=value;key2=value;...".
    infos = fields[7].split(';')

    for i, info in enumerate(infos, 1):
        # info should be "key=value".
        try:
            key, value = info.split('=')
        # But sometimes it is just "value", so we'll make our own key.
        except ValueError:
            # key = 'INFO{}'.format(i)
            # value = info
            continue
        # Set the value to None if there is no value.
        result[key] = _get_value(value)
    # print(result)
    return result


def _get_value(value):
    """Interpret null values and return ``None``. Return a list if the value
    contains a comma.
    """
    if not value or value in ['', '.', 'NA']:
        return None
    if ',' in value:
        return value.split(',')
    return value


def _count_comments(filename):
    """Count comment lines (those that start with "#") in an optionally
    gzipped file.
    :param filename:  An optionally gzipped file.
    """
    comments = 0
    fn_open = gzip.open if filename.endswith('.gz') else open
    with fn_open(filename) as fh:
        for line in fh:
            if line.startswith('#' or '##'):
                comments += 1
            else:
                break
    return comments


def _vcf_filter(df, *argv):
    arguement = " & ".join(argv)
    data = df.query(arguement)
    return data


data = dataframe(filename, large=False)
final_data = add_samples(filename, data)

# This statement return the location of the FORMAT column in vcf file
format_index = (final_data.columns.get_loc("FORMAT"))

# This statement return the location of index of the first sample
first_Sample_index = format_index + 1

# This statement returns the length of the data frame
Num_cols = len(final_data.columns)

# Thi statement returns the sample names
sample_names = final_data.columns.values

# This method is parsing the sample field and returning only genotype field
for i in range(first_Sample_index, Num_cols):
    s = sample_names[i]
    final_data[s] = final_data[s].str.split(':').str[0]



final_data.to_csv(output_file)
