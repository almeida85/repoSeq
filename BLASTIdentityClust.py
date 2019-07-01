#!/usr/bin/env python
#
#       Identity-based Sequences Clustering.
#       Perform a BLAST pairwise alignments between two groups of sequences and
#       group them according to an identity value.
#
#       Copyright 2010 Yasser Almeida Hernandez <almeida@cim.sld.cu>
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.
#
#
#

import sys, os, time
from commands import getstatusoutput

__doc__='''
        Identity-based Sequences Clustering.
        Perform a BLAST pairwise alignments between two groups of sequences and
        group them according to an identity value.
        '''

log_file = open('seqs_clustering.log', 'a+')

class IdentityCutoffError(Exception):
    """Identity cutoff out of range..."""
    pass

class Fasta:
    '''
    Fasta class. Define sequences in FASTA format
    '''
    def __init__(self, fasta_name):
        '''
        Fasta Atributtes

        self.fasta_dict = Dictionary to store the .fasta file, where the key is
                          the sequence's name and the value is the sequence.
        self.seqnames_list = List of sequences names in the .fasta file (keys).
        self.nseqs = Number of sequences.
        log_file = Verbose log file.
        '''

        self.fasta_dict = {}
        self.fasta_dict = self.read_fasta(fasta_name)
        self.seqnames_list = []
        self.seqnames_list = self.fasta_dict.keys()
        self.nseqs = len(self.seqnames_list)
        #log_file = open('seqs_clustering.log', 'a+')

    def read_fasta(self, fasta_name):
        '''
        Read  and store the .fasta file in a dictionary.
        '''

        fasta_file = open(fasta_name, "r")
        sequence_name=''
        cant_seqs = 0
        seq_dict = {}

        for line in fasta_file:
            line=line.strip()
            if line:
                if line[0] == '>':
                    sequence_name=line[1:]
                    seq_dict[sequence_name] = ''
                    cant_seqs += 1
                else:
                    if not sequence_name:
                        print 'Ignoring line',line
                    else:
                        seq_dict[sequence_name] += line

        print cant_seqs,'sequences read from',fasta_name
        log_file.write(str(cant_seqs) + ' sequences read from ' + fasta_name + '\n\n')
        fasta_file.close()
        return seq_dict

    def blast(self, fasta_name, reference_fasta_name):
        '''
        Run BLAST.
        '''

        ## Formating the reference sequences database file.
        print '\nFormating the reference sequences database file...'
        log_file.write('Formating the reference sequences database file...\n')
        formatdb_string = 'formatdb -t REFERENCE -i "%s"' % reference_fasta_name
        os.system(formatdb_string)

        ## Running BLAST.
        print "\nRunning BLAST..."
        log_file.write("\nRunning BLAST...\n")
        blast_cmd = 'blast2 -d "%s" -i "%s" -p blastp -m 8 | cut -f1,2,3 > blast.res' % (reference_fasta_name, fasta_name)
        st, out = getstatusoutput(blast_cmd)
        print "BLAST DONE\n\n"
        log_file.write("BLAST DONE\n\n")

        return 0

    def read_identity_file(self,identity_file):
        '''
        Read the BLAST output results file.
        '''

        print "Reading BLAST results..."
        log_file.write("Reading BLAST results...\n")

        identity_list = []               # List of lists with the query and subject sequences names
        query_seq_name = ''              # Query ligand name
        subject_seq_name = ''            # Reference ligand name

        for line in identity_file:
            line=line.strip()
            line_words = line.split()
            query_seq_name = line_words[0]
            subject_seq_name = line_words[1]
            identity_value = float(line_words[2])
            identity_list.append((query_seq_name,subject_seq_name,identity_value))

        identity_file.close()
        return identity_list

    def identity_clustering(self, identity_list, reference_fasta_name, cutoff,):
        '''
        Perform the sequences identity clustering.
        First, eliminate the redundance and them group the sequences.
        '''

        ## Reference sequences file.
        RefFasta = Fasta(reference_fasta_name)

        ## Filtering sequences aligned with itself...
        sequences_aligned = []           # Names of aligned sequences.

        for items in identity_list:
            sequences_aligned.append(items[0])

        print "-> Ignoring sequences with any alignment..."
        log_file.write("Ignoring sequences with any alignment...\n")

        for name in self.seqnames_list:
            if sequences_aligned.count(name) == 1:
                self.seqnames_list.remove(name)
                log_file.write('-> ' + name + ' only produce alignments with itself.\n')
                for items in identity_list:
                    if items[0].rsplit('|', 1)[:1] == name and items[1].rsplit('|', 1)[:1] == name:
                        identity_list.remove(items)

        ## Ignoring redundance...
        print "-> Ignoring redundance..."
        log_file.write("\nIgnoring redundance...\n")
        print '-> Writing sequences clusters with identity >= of',cutoff,'%...\n'
        log_file.write('Writing sequences clusters with identity >= of ' + str(cutoff) + '%...\n')

        ignored_seqs_list = []           # Sequences ignored (identity = 100%) for eliminate redundance

        for seq in self.seqnames_list:
            if seq in ignored_seqs_list:
                continue
            else:
                fasta_cluster = open(seq + "_cluster.fasta", "w")
                fasta_cluster.write('>' + seq + '|CLUSTER_HEAD\n' + self.fasta_dict[seq] + '\n')
                log_file.write('\n\nCluster head: ' + seq + '\n')

                for item in identity_list:
                    if item[0] == seq:
                        identity = item[2]

                        if item[0] == item[1] and identity == 100:
                            continue
                        elif item[0] != item[1] and identity == 100:
                            ignored_seqs_list.append(item[1])
                        elif identity < cutoff:
                            continue
                        elif item[0] != item[1] and cutoff <= identity <= 100:
                            fasta_cluster.write('>' + item[1] + '|' + str(round(identity, 2)) + '\n' + RefFasta.fasta_dict[item[1]] + '\n')
                            log_file.write(item[1] + '\t')
                    else:
                        continue

        return 0

def print_duration(seconds):
    output=''
    seconds_per_minute=60
    seconds_per_hour=60*seconds_per_minute
    seconds_per_day=24*seconds_per_hour

    days=int(seconds/seconds_per_day)
    hours=int((seconds % seconds_per_day)/seconds_per_hour)
    minutes=int((seconds % seconds_per_hour)/seconds_per_minute)
    seconds=seconds%60

    if days:
        output+=('%d' % (days))
        if days==1:
            output+=' day '

        else:
            output+=' days '

    if hours:
        output+=('%2d' % (hours))
        if hours==1:
            output+=' hour '

        else:
            output+=' hours '

    if minutes:
        output+=('%2d' % (minutes))
        if minutes==1:
            output+=' minute '

        else:
            output+=' minutes '

    if seconds:
        output+=('%4.2f' % (seconds))
        if seconds==1:
            output+=' second '

        else:
            output+=' seconds '

    return output

def usage():
    '''
    HELP
    '''

    print 'USAGE: %s <query_fasta_file> <identity_cut-off>\n' % sys.argv[0]
    log_file.write('\nERROR!!: Argument missing...\n')
    print ''
    return 1

def run():
    '''
    RUN THE SCRIPT
    '''

    if len(sys.argv)!= 3:
        print '\nERROR!!: Argument missing...'
        usage()
        sys.exit()
    else:
        fasta_name = sys.argv[1]
        #reference_fasta_name = sys.argv[2]
        
        ## PDB 'apo' sequences...
        reference_fasta_name = "/home/almeida/Travail/Residues_fiting/fastas/seqs_without_ligands/seqs_without_ligands_PDB_query_contacts_purged.fasta"
        
        identity_cutoff = float(sys.argv[2])

        if identity_cutoff < 0 or identity_cutoff >= 100:
            error_test='\n\nIdentity cut-off out of range!!\n'
            log_file.write('ERROR!!: Identity cut-off out of range.\n\n')

            print '\n', time.ctime(time.time())
            log_file.write(time.ctime(time.time()))

            raise IdentityCutoffError(error_test)


        ## The object Fasta take at least two parameters.
        ## 'fasta_name' is the main sequence file name.'
        ## 'reference_fasta_name' is the subject sequence file name database for
        ## run the BLAST program. If 'fasta_name' and 'reference_fasta_name' are
        ## the same (i mean, if you will run the BLAST against the main sequence
        ## file it self), Fasta take only one parameters: 'fasta_name'.
        MFasta = Fasta(fasta_name)

        ## Call BLAST
        MFasta.blast(fasta_name, reference_fasta_name)

        ## Read the BLAST output and perform the clustering
        identity_file = open("blast.res", "r")
        identity_list = MFasta.read_identity_file(identity_file)
        MFasta.identity_clustering(identity_list, reference_fasta_name, identity_cutoff)

        return 0

# MAIN
if __name__ == "__main__":
    log_file = open('seqs_clustering.log', 'w')

    print '*************************************'
    print '*** SEQUENCES IDENTITY CLUSTERING ***'
    print '*************************************'

    log_file.write('*************************************\n')
    log_file.write('*** SEQUENCES IDENTITY CLUSTERING ***\n')
    log_file.write('*************************************\n')

    print time.ctime(time.time())
    initial_time = time.time()
    log_file.write(time.ctime(time.time()) + '\n\n')

    run()

    final_time=time.time()

    print 'Script finished in %s' % print_duration(final_time-initial_time)
    log_file.write('\nScript finished in ' + print_duration(final_time-initial_time) + '\n')

    print 'DONE'
    log_file.write('DONE\n')

    print time.ctime(time.time()),'\n'
    log_file.close()
