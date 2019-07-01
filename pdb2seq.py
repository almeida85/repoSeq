#!/usr/bin/python

"""
Extract the sequence of a .pdb file and optionally, run a local
BLAST search with the sequence(s).
You need the BioPython module to run this script

Author: Yasser Almeida Hdez.
yasser.almeida@gmail.com
Version 1.0
August 2016
"""

import string, sys, optparse
from Bio.PDB import *
from Bio.SeqIO import *

## Classes
class BioPDBMissing(Exception):
    """Bio.PDB Module missing..."""
    pass

class PDBNotFound(Exception):
    """PDB file not found..."""
    pass

class SequenceBuilder:
    def __init__(self, pdb_id, pdb_filename, query_chains=None):
        self.structure_id = pdb_id
        self.structure = self._makeStructure(pdb_filename)
        self.structure_chains = self._obtain_chains()
        self.chains_to_sequence = self._check_chains(query_chains)


    ## Private methods
    def _makeStructure(self, pdb_file_name):
        '''Make the stucture objet...'''

        structure = PDBParser()
        structure = structure.get_structure('crystal', pdb_file_name)
        return structure

    def _obtain_chains(self):
        '''Get the structure chains'''

        chain_list = []
        for chain in self.structure.get_chains():
            if chain.get_id() == " ":
                continue
            else:
                chain_list.append(chain.get_id())

        return chain_list

    def _check_chains(self, query_chains):
        '''Check if the queries chains exist in the structure and obtain the sequences'''

        if query_chains == None:
           query_chains = self.structure_chains
           return query_chains
        else:
            chains_to_sequence = []
            query_chains = list(query_chains)

            for chain in query_chains:
                if chain not in self.structure_chains:
                    print "Chain",chain,"do not exist in the structure."
                    continue
                else:
                    chains_to_sequence.append(chain)

            return chains_to_sequence


    ## Public methods
    def get_sequence(self, redundancy):
        '''Make the FASTA file from the PDB file...'''

        chains_sequences = {}

        ## Iterating over protein chains...
        for chain in self.structure.get_chains():
            chain_id = chain.get_id()
            if chain_id not in self.chains_to_sequence:
                continue
            elif chain_id == " ":
                continue
            else:
                res_list = []
                ## Iterating over residues in chains...
                for residue in chain.get_list():
                    res_id = residue.get_resname()
                    ## Translating to 1 letter code...
                    if to_one_letter_code.has_key(res_id) == True:
                        r = to_one_letter_code.get(res_id)
                        res_list.append(r)
                    else:
                        continue

                sequence = string.join(res_list, sep='')

                if redundancy == True:
                    if sequence in chains_sequences.values():
                        continue
                    else:
                        chains_sequences[chain_id] = ''
                        chains_sequences[chain_id] = sequence
                else:
                    chains_sequences[chain_id] = ''
                    chains_sequences[chain_id] = sequence

        return chains_sequences

    def write_sequence(self, sequences, format, unique_sequence):
        '''Format the sequence...'''

        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.Alphabet import IUPAC

        if unique_sequence == True:
            for chain in sequences.keys():
                output_filename = self.structure_id + "_" + chain + "." + format
                output_file = open(output_filename, "w")

                record = SeqRecord(Seq(sequences[chain],IUPAC.protein), id=self.structure_id + '|Chain_' + chain, name=chain,\
                                   description=self.structure.header['name'])

                output_file.write(record.format(format) + '\n')
                output_file.close()

            return 0
        else:
            if self.chains_to_sequence != self.structure_chains:
                output_filename = self.structure_id + "_" + string.join(self.chains_to_sequence, sep="") +"." + format
            else:
                output_filename = self.structure_id + "." + format

            output_file = open(output_filename, "w")

            for chain in sequences.keys():
                record = SeqRecord(Seq(sequences[chain],IUPAC.protein), id=self.structure_id + '|Chain_' + chain, name=chain,\
                                   description=self.structure.header['name'])

                output_file.write(record.format(format) + '\n')

            output_file.close()

            return output_filename

    def blast(self, database, input_file, unique_sequence):
        '''
        Run local BLAST.
        '''

        import os
        from commands import getstatusoutput

        if unique_sequence == True:
            parse.error('Unset the -u option to run BLAST')
            return -1

        ## Formating the reference sequences database file.
        print 'Formating the sequences database file...'
        formatdb_string = 'formatdb -t REFERENCE -i "%s"' % database
        os.system(formatdb_string)

        ## Running BLAST.
        print "Running BLAST..."
        blast_cmd = 'blast2 -d "%s" -i "%s" -p blastp -m 9 > "%s"_blast.res' % \
                       (database, input_file, self.structure_id)
        st, out = getstatusoutput(blast_cmd)
        print "BLAST DONE"

        return

    def stats(self, sequences):
        '''Generate a file with useful data of the protein. Based in ProtParam Tools from Expasy'''

        from Bio.SeqUtils import ProtParam

        #output_filename = self.structure_id + "_stats.dat"

        #stats_file = open(output_filename, "w")

        for chain in sequences.keys():
            print "*************************"
            print "**  ", self.structure_id,"- Chain",chain,"   **"
            print "*************************"

            seq_stats = ProtParam.ProteinAnalysis(sequences[chain])


            ## Printing the aa count...
            print "Amino acids counts and percents"
            total_aa = 0
            for aa, percent in zip(seq_stats.count_amino_acids(), seq_stats.get_amino_acids_percent()):
                print aa,":",seq_stats.count_amino_acids()[aa],"(",round(seq_stats.get_amino_acids_percent()[aa]*100, 2),"% )"
                total_aa += seq_stats.count_amino_acids()[aa]
            print "TOTAL:", total_aa

            molar_mass = seq_stats.molecular_weight()/100

            print "\nMolecular mass:",molar_mass,"kDa"

            print "\nIsoelectric point:", round(seq_stats.isoelectric_point(), 2)

            extintion_coef_Cyst, extintion_coef_noCyst, molar_extintion_coef_Cyst, molar_extintion_coef_noCyst = self.get_extintion_coef(seq_stats, molar_mass)

            print "\nExtintion coefficient (Cystines) =", extintion_coef_Cyst,"M^-1*cm^-1"
            print "Extintion coefficient (no Cystines) =", extintion_coef_noCyst,"M^-1*cm^-1"

            print "\nMolar extintion coefficient (Cystines) [Abs 0.1% (=1 g/l)] =", molar_extintion_coef_Cyst
            print "Molar extintion coefficient (no Cystines) [Abs 0.1% (=1 g/l)] =", molar_extintion_coef_noCyst,"\n"

        return

    def get_extintion_coef(self, sequence, molar_mass):

        ext_Tyr = 1490
        ext_Trp = 5500
        ext_Cyst = 125

        nTyr = sequence.count_amino_acids()['Y']
        nTrp = sequence.count_amino_acids()['W']
        nCys = sequence.count_amino_acids()['C']
        nCyst = nCys/2

        extintion_coef_Cyst = nTyr*ext_Tyr + nTrp*ext_Trp + nCyst*ext_Cyst
        extintion_coef_noCyst = nTyr*ext_Tyr + nTrp*ext_Trp

        molar_extintion_coef_Cyst = round(extintion_coef_Cyst/molar_mass, 2) / 100
        molar_extintion_coef_noCyst = round(extintion_coef_noCyst/molar_mass, 2) / 100

        return extintion_coef_Cyst, extintion_coef_noCyst, molar_extintion_coef_Cyst, molar_extintion_coef_noCyst



## Functions
def checkPDBFile(pdb_filename):
    '''Checking the input pdb file...'''
    try:
        infile = open(pdb_filename, "r")
        infile.close()
    except:
        usage()
        error_text = 'PDB file not found!!!'
        raise PDBNotFound(error_text)

def run_script():
    '''
    This function run the script...
    '''

    usage = "%prog [options] PDB_file"
    description = "This program get the sequence(s) of protein chain(s) from a PDB file."
    version = "1.0"

    parser = optparse.OptionParser(usage=usage,description=description,version=version)

    parser.add_option("-f", dest = "format", type="choice", choices = ['fasta', 'embl', 'gb', 'tab', 'imgt'], default = "fasta", \
                      help = "Output sequence format: fasta, embl, gb, imgt, tab [default: %default].")
    parser.add_option("-c", dest = "chains", help = "Specify chain(s) to extract sequence. Type the chains ids with no space (ej. -c AB).")
    parser.add_option("-r", dest = "remove_redundancy", action = "store_true", default = False, help = "Remove redundant sequences (ej. homodimers, trimers, etc.)")
    parser.add_option("-u", dest = "unique_sequence", action = "store_true", default = False, help = "Write each chain sequence in a file independently.")
    parser.add_option('-d', dest = "database", help = "Sequence database for local BLAST search.")
    parser.add_option('--blast', dest = "blast", action = "store_true", default = False, help = "Run a local BLAST search.")
    parser.add_option('--stats', dest = "stats", action = "store_true", help = "Calculate diferent properties of the protein.")

    (options, args) = parser.parse_args()

    if not len(args):
        parser.error('You need to provide a .pdb file.\nSee help (-h, --help) for more options.')

    if options.blast and options.database == None:
        parser.error('You must to set the sequence database file to run BLAST.')

    pdb_filename = sys.argv[1]
    pdb_id = pdb_filename.rsplit('.',1)[0]

    ## Check if the .pdb file exist...
    checkPDBFile(pdb_filename)

    PDBsequence = SequenceBuilder(pdb_id, pdb_filename, options.chains)

    sequences = PDBsequence.get_sequence(options.remove_redundancy)

    #print sequences

    seq_filename = PDBsequence.write_sequence(sequences, options.format, options.unique_sequence)

    if options.blast:
        PDBsequence.blast(options.database, seq_filename, options.unique_sequence)

    if options.stats == True:
        PDBsequence.stats(sequences)

    return 0

# Main
if __name__ == "__main__":
    run_script()
