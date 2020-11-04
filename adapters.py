"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Porechop

This module contains the class and sequences for known adapters used in Oxford Nanopore library
preparation kits.

This file is part of Porechop. Porechop is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Porechop is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Porechop. If
not, see <http://www.gnu.org/licenses/>.


##### MODIFIED VERSION FOR EHD SAMPLES ########## 

JDGJ

#I treat p5 and p7 as adapters to be trimmed.
#I treat FLDXXXX and CS1/CS2 as barcodes for sorting.
#MinION sequences in both directions, so include end_sequence as reverse compliment for all.
#I also commented out out the 1D^2 adapters since they were being trimmed from some reads even though
I didn't use the kits for library prep. SQK-NSK007 adapters are the correct adapters to be trimmed
for the kit I used.
#To look at pX/csX associations, I did two runs of Porechop where I sorted reads into p5 and p7 first,
then a second run to sort reads in each one into cs1 and cs2 folders. I then looked at number of reads 
in each directory /pX/csX
#For demultiplexing, I removed the exisiting ONT barcodes in the original script to avoid conflicts,
then I added FLDXXXX barcodes
"""


class Adapter(object):

    def __init__(self, name, start_sequence=None, end_sequence=None, both_ends_sequence=None):
        self.name = name
        self.start_sequence = start_sequence if start_sequence else []
        self.end_sequence = end_sequence if end_sequence else []
        if both_ends_sequence:
            self.start_sequence = both_ends_sequence
            self.end_sequence = both_ends_sequence
        self.best_start_score, self.best_end_score = 0.0, 0.0

    def best_start_or_end_score(self):
        return max(self.best_start_score, self.best_end_score)

    def is_barcode(self):
        return self.name.startswith('Barcode ')

    def barcode_direction(self):
        if '_rev' in self.start_sequence[0]:
            return 'reverse'
        else:
            return 'forward'
    

    def get_barcode_name(self):
        """
        Gets the barcode name for the output files. We want a concise name, so it looks at all
        options and chooses the shortest.
        """
        possible_names = [self.name]
        if self.start_sequence:
            possible_names.append(self.start_sequence[0])
        if self.end_sequence:
            possible_names.append(self.end_sequence[0])
        barcode_name = sorted(possible_names, key=lambda x: len(x))[0]
        return barcode_name.replace(' ', '_')


# INSTRUCTIONS FOR ADDING CUSTOM ADAPTERS
# ---------------------------------------
# If you need Porechop to remove adapters that aren't included, you can add your own my modifying
# the ADAPTERS list below.
#
# Here is the format for a normal adapter:
#     Adapter('Adapter_set_name',
#             start_sequence=('Start_adapter_name', 'AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT'),
#             end_sequence=('End_adapter_name', 'AACCGGTTAACCGGTTAACCGGTTAACCGGTT'))
#
# You can exclude start_sequence and end_sequence as appropriate.
#
# If you have custom Barcodes, make sure that the adapter set name starts with 'Barcode '. Also,
# remove the existing barcode sequences from this file to avoid conflicts:
#     Adapter('Barcode 1',
#             start_sequence=('Barcode_1_start', 'AAAAAAAACCCCCCCCGGGGGGGGTTTTTTTT'),
#             end_sequence=('Barcode_1_end', 'AAAAAAAACCCCCCCCGGGGGGGGTTTTTTTT')),
#     Adapter('Barcode 2',
#             start_sequence=('Barcode_2_start', 'TTTTTTTTGGGGGGGGCCCCCCCCAAAAAAAA'),
#             end_sequence=('Barcode_2_end', 'TTTTTTTTGGGGGGGGCCCCCCCCAAAAAAAA'))


ADAPTERS = [Adapter('Barcode 105 (forward)',
                    start_sequence=('Barcode_105_start', ''),
                    end_sequence=('Barcode_105_end', 'TGATAGAGAG')),
            Adapter('Barcode 106 (forward)',
                    start_sequence=('Barcode_106_start', ''),
                    end_sequence=('Barcode_106_end', 'GCTACTAGCG')),
            Adapter('Barcode 107 (forward)',
                    start_sequence=('Barcode_107_start', ''),
                    end_sequence=('Barcode_107_end', 'TGCGAGACGT')),
            Adapter('Barcode 108 (forward)',
                    start_sequence=('Barcode_108_start', ''),
                    end_sequence=('Barcode_108_end', 'CGATGACAGA')),
            Adapter('Barcode 109 (forward)',
                    start_sequence=('Barcode_109_start', ''),
                    end_sequence=('Barcode_109_end', 'GACTCATGCT')),
            Adapter('Barcode 110 (forward)',
                    start_sequence=('Barcode_110_start', ''),
                    end_sequence=('Barcode_110_end', 'GTCTGATACG')),
            Adapter('Barcode 111 (forward)',
                    start_sequence=('Barcode_111_start', ''),
                    end_sequence=('Barcode_111_end', 'ACTAGCTGTC')),
            Adapter('Barcode 112 (forward)',
                    start_sequence=('Barcode_112_start', ''),
                    end_sequence=('Barcode_112_end', 'GCGTAGACGA')),
            Adapter('Barcode 113 (forward)',
                    start_sequence=('Barcode_113_start', ''),
                    end_sequence=('Barcode_113_end', 'CTCAGCAGTG')),
            Adapter('Barcode 114 (forward)',
                    start_sequence=('Barcode_114_start', ''),
                    end_sequence=('Barcode_114_end', 'CAGTCTACAT')),
            Adapter('Barcode 115 (forward)',
                    start_sequence=('Barcode_115_start', ''),
                    end_sequence=('Barcode_115_end', 'TACTGCAGCG')),
            Adapter('Barcode 116 (forward)',
                    start_sequence=('Barcode_116_start', ''),
                    end_sequence=('Barcode_116_end', 'TACACAGTAG')),
            Adapter('Barcode 117 (forward)',
                    start_sequence=('Barcode_117_start', ''),
                    end_sequence=('Barcode_117_end', 'CACATACAGT')),
            Adapter('Barcode 118 (forward)',
                    start_sequence=('Barcode_118_start', ''),
                    end_sequence=('Barcode_118_end', 'CACAGTGATG')),
            Adapter('Barcode 119 (forward)',
                    start_sequence=('Barcode_119_start', ''),
                    end_sequence=('Barcode_119_end', 'CGAGCTAGCA')),
            Adapter('Barcode 120 (forward)',
                    start_sequence=('Barcode_120_start', ''),
                    end_sequence=('Barcode_120_end', 'GAGACTATGC'))]
           
            


def make_full_native_barcode_adapter(barcode_num):
    barcode = [x for x in ADAPTERS if x.name == 'Barcode ' + str(barcode_num) + ' (reverse)'][0]
    start_barcode_seq = barcode.start_sequence[1]
    end_barcode_seq = barcode.end_sequence[1]

    start_full_seq = 'AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAA' + start_barcode_seq + 'CAGCACCT'
    end_full_seq = 'AGGTGCTG' + end_barcode_seq + 'TTAACCTTAGCAATACGTAACTGAACGAAGT'

    return Adapter('Native barcoding ' + str(barcode_num) + ' (full sequence)',
                   start_sequence=('NB' + '%02d' % barcode_num + '_start', start_full_seq),
                   end_sequence=('NB' + '%02d' % barcode_num + '_end', end_full_seq))


def make_old_full_rapid_barcode_adapter(barcode_num):  # applies to SQK-RBK001
    barcode = [x for x in ADAPTERS if x.name == 'Barcode ' + str(barcode_num) + ' (forward)'][0]
    start_barcode_seq = barcode.start_sequence[1]

    start_full_seq = 'AATGTACTTCGTTCAGTTACG' + 'TATTGCT' + start_barcode_seq + \
                     'GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA'

    return Adapter('Rapid barcoding ' + str(barcode_num) + ' (full sequence, old)',
                   start_sequence=('RB' + '%02d' % barcode_num + '_full', start_full_seq))


def make_new_full_rapid_barcode_adapter(barcode_num):  # applies to SQK-RBK004
    barcode = [x for x in ADAPTERS if x.name == 'Barcode ' + str(barcode_num) + ' (forward)'][0]
    start_barcode_seq = barcode.start_sequence[1]

    start_full_seq = 'AATGTACTTCGTTCAGTTACG' + 'GCTTGGGTGTTTAACC' + start_barcode_seq + \
                     'GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA'

    return Adapter('Rapid barcoding ' + str(barcode_num) + ' (full sequence, new)',
                   start_sequence=('RB' + '%02d' % barcode_num + '_full', start_full_seq))
