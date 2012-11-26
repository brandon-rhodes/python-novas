#!/usr/bin/env python

from __future__ import with_statement

import os
import re
import struct

def download_ascii(ascii_dir, de_number, verbose=True):
    import ftplib
    
    files = []
    
    if verbose:
        print("Retrieving list of ephemeris files")

    jpl = ftplib.FTP('ssd.jpl.nasa.gov')
    jpl.login()
    jpl.cwd('pub/eph/planets/ascii/de%s' % de_number)
    jpl.retrlines('LIST', files.append)
    
    for filename in [file.split()[-1] for file in files]:
        filelines = []
        if verbose:
            print("Downloading ASCII ephemeris file %s" % filename)
        jpl.retrlines('RETR %s' % filename, filelines.append)
        outfile = open(os.path.join(ascii_dir, filename), 'w')
        for line in filelines:
            outfile.write(line + '\n')
        outfile.close()
        
    jpl.close()

def process_header(ascii_dir, de_number, binary_file, verbose=True):
    group_nums = [0, 1010, 1030, 1040, 1041, 1050, 1070]
    
    header_file = open(os.path.join(ascii_dir, "header.%s" % de_number), 'r')
    header_lines = header_file.read()
    header_file.close()

    header = [group.strip(' ' + os.linesep) for group in re.split(r'GROUP\s*\d{4}', header_lines)]

    groups = dict(zip(group_nums, header))

    ksize, ncoeff = [int(number) for number in re.findall(r'\d+', groups[0])]

    alphanumeric_records = [line.strip() for line in groups[1010].split(os.linesep)]

    start = float(groups[1030].split()[0])
    end = float(groups[1030].split()[1])
    span = float(groups[1030].split()[2])

    ncon = int(groups[1040].split()[0])
    constants = dict(zip(groups[1040].split()[1:], [float(constant.replace('D', 'e')) for constant in groups[1041].split()[1:]]))

    au = float(constants['AU'])
    emrat = float(constants['EMRAT'])
    denum = int(constants['DENUM'])

    pointers = [int(pointer) for pointer in groups[1050].split()]

    for line in alphanumeric_records:
        binary_file.write(line.ljust(84))

    for key in constants.keys():
        binary_file.write(key.ljust(6))

    padding = 2400 - ncon * 6

    binary_file.write(struct.pack("=%ix" % padding))

    binary_file.write(struct.pack('=dddidd', start, end, span, ncon, au, emrat))

    for i in range(12):
        binary_file.write(struct.pack('=iii', pointers[i], pointers[i+13], pointers[i+13+13]))

    binary_file.write(struct.pack('=i', denum))

    binary_file.write(struct.pack('=iii', pointers[12], pointers[12+13], pointers[12+13+13]))

    binary_file.write(struct.pack('=5288x'))

    for key in constants.keys():
        binary_file.write(struct.pack('=d', constants[key]))

    padding = 8144 - ncon * 8

    binary_file.write(struct.pack('=%ix' % padding))
    
    return ncoeff
    
def process_data_files(data_files, ncoeff, binary_file, verbose=True):
    old_block = []
    new_block = []
    
    for file in data_files:
        if verbose:
            print("Processing ASCII ephemeris file %s" % os.path.split(file)[-1])
        with open(file, 'r') as data_file:
            for line in data_file:
                if len(line.split()) == 2:
                    if new_block != old_block:
                        for number in new_block[:ncoeff]:
                            binary_file.write(struct.pack('=d', number))
                        old_block = new_block
                        new_block = []
                    else:
                        new_block = []
                else:
                    for number in line.replace('D', 'e').split():
                        new_block.append(float(number))

    for number in new_block[:ncoeff]:
        binary_file.write(struct.pack('=d', number))
        
if __name__ == '__main__':
    import tempfile
    import shutil
    
    de_number = 405
    
    current_dir = os.getcwd()
    temp_dir = tempfile.mkdtemp()
    
    try:
        download_ascii(temp_dir, de_number)
        
        binary_file = open(os.path.join(current_dir, "DE%s.bin") % de_number, 'wb')
        
        ncoeff = process_header(temp_dir, de_number, binary_file)
        
        data_files = [os.path.join(temp_dir, "ascp%d.%s" % (year, de_number)) for year in range(1600, 2220, 20)]
        
        process_data_files(data_files, ncoeff, binary_file)
        
        binary_file.close()
    finally:
        shutil.rmtree(temp_dir)
