import os
import sqlite3
from bioagents import BioagentException
import csv

tcga_study_names = ['ACC', 'BLCA', 'BRCA', 'CESC','CHOL', 'COAD', 'COADREAD', 'DLBC', 'GBM', 'GBMLGG', 'HNSC',
                    'KICH', 'KIPAN','KIRC', 'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'OV', 'PAAD', 'PCPG',
                    'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'STES', 'TGCT', 'THCA', 'UCEC', 'UCS', 'UVM']
class DatabaseInitializer:
    """ Fills the pnnl database from the given data files"""

    def __init__(self, path):
        db_file = os.path.join(path, 'pnnl-dataset.db')
        if os.path.isfile(db_file):
            self.cadb = sqlite3.connect(db_file)
        else:
            # create table if it doesn't exist
            fp = open(db_file, 'w')
            fp.close()
            self.cadb = sqlite3.connect(db_file)
            self.populate_tables(path)

    def __del__(self):
        self.cadb.close()

    def populate_tables(self, path):
        """
        Fills all the tables in the database
        :param path: Path to the folder that keeps all the data files
        :return:
        """
        self.populate_correlation_table(path)
        self.populate_causality_table(path)
        self.populate_mutsig_table(path)
        self.populate_unexplained_table()
        self.populate_explained_table()
        self.populate_sif_relations_table(path)
        self.populate_mutex_table(path)
        self.populate_tcga_names_table(path)

    def populate_causality_table(self, path):
        """
        Fills the causality table
        :param path: Path to the folder that keeps causative-data-centric.sif
        :return:
        """
        opposite_rel = {
            'phosphorylates': 'is-phosphorylated-by',
            'dephosphorylates': 'is-dephosphorylated-by',
            'upregulates-expression': 'expression-is-upregulated-by',
            'downregulates-expression': 'expression-is-downregulated-by',
        }

        try:
            causality_path = os.path.join(path, 'causative-data-centric.sif')
        except Exception as e:
            raise BioagentException.PathNotFoundException()

        causality_file = open(causality_path, 'r')

        with self.cadb:
            cur = self.cadb.cursor()
            cur.execute("DROP TABLE IF EXISTS Causality")
            cur.execute("CREATE TABLE Causality(Id1 TEXT, PSite1 TEXT, Id2 TEXT, PSite2 TEXT, Rel TEXT, UriStr TEXT)")

            for line in causality_file:
                vals = line.split('\t')
                id_str1 = vals[0].upper().split('-')
                id1 = id_str1[0]
                if len(id_str1) > 1:
                    p_site1 = id_str1[1]
                else:
                    p_site1 = ' '

                id_str2 = vals[2].upper().split('-')
                id2 = id_str2[0]

                if len(id_str2) > 1:
                    p_site2 = id_str2[1]
                else:
                    p_site2 = ' '

                rel = vals[1]

                uri_arr = []
                if vals[3]:
                    uri_arr = vals[3].split(" ")

                if len(uri_arr) == 0:
                    uri_arr = [vals[3]]

                uri_str = ""
                for uri in uri_arr:
                    uri_str = uri_str + "uri= " + uri + "&"

                opp_rel = opposite_rel[rel]
                cur.execute("INSERT INTO Causality VALUES(?, ?, ?, ?, ?, ?)", (id1, p_site1, id2, p_site2, rel, uri_str))
                # opposite relation
                cur.execute("INSERT INTO Causality VALUES(?, ?, ?, ?, ?, ?)", (id2, p_site2, id1, p_site1, opp_rel, uri_str))

        causality_file.close()

    def populate_correlation_table(self, path):
        """
        Fills the correlation table
        :param::path: Path to the folder that keeps PNNL-ovarian-correlations.txt
        :return:
        """

        try:
            pnnl_path = os.path.join(path, 'PNNL-ovarian-correlations.txt')
        except Exception as e:
            raise BioagentException.PathNotFoundException()
        pnnl_file = open(pnnl_path, 'r')

        with self.cadb:
            cur = self.cadb.cursor()
            cur.execute("DROP TABLE IF EXISTS Correlations")
            cur.execute("CREATE TABLE Correlations(Id1 TEXT, PSite1 TEXT, Id2 TEXT, PSite2 TEXT, Corr REAL, PVal REAL)")

            for line in pnnl_file:
                if line.find('/') > -1:  # incorrectly formatted strings
                    continue
                vals = line.split('\t')
                id_str1 = vals[0].upper().split('-')
                id1 = id_str1[0]
                if len(id_str1) > 1:
                    p_site1 = id_str1[1]
                else:
                    p_site1 = ' '

                id_str2 = vals[1].upper().split('-')
                id2 = id_str2[0]

                if len(id_str2) > 1:
                    p_site2 = id_str2[1]
                else:
                    p_site2 = ' '

                corr = float(vals[2].rstrip('\n'))

                p_val = float(vals[3].rstrip('\n'))

                cur.execute("INSERT INTO Correlations VALUES(?, ?, ?, ?, ?, ?)", (id1, p_site1, id2, p_site2, corr, p_val))

        pnnl_file.close()

    def populate_mutsig_table(self, path):
        """
        :param path: Path to the folder that keeps TCGA folder and mutsig.txt files
        :return:
        """

        try:
            mutsig_path = os.path.join(path, 'TCGA')
        except Exception as e:
            raise BioagentException.PathNotFoundException()

        folders = os.listdir(mutsig_path)

        with self.cadb:
            cur = self.cadb.cursor()
            cur.execute("DROP TABLE IF EXISTS MutSig")
            cur.execute("CREATE TABLE MutSig(Id TEXT, Disease TEXT, PVal REAL, QVal Real)")

            for folder in folders:
                try:
                    if folder in tcga_study_names:
                        disease_path = os.path.join(mutsig_path, folder)
                        file_path = os.path.join(disease_path, 'scores-mutsig.txt')
                    else:
                        continue
                except Exception as e:
                    raise BioagentException.PathNotFoundException()

                mutsig_file = open(file_path, 'r')
                next(mutsig_file) # skip the header line

                for line in mutsig_file:
                    vals = line.split('\t')
                    gene_id = vals[1]
                    p_val = vals[17]
                    q_val = vals[18].rstrip('\n')
                    cur.execute("INSERT INTO MutSig VALUES(?, ?, ?, ?)", (gene_id, folder,  p_val, q_val))

                mutsig_file.close()

    def populate_mutex_table(self, path):
        """
        Finds mutually exclusive gene groups
        :param path: Path to the folder that keeps ranked-groups.txt
        :return:
        """

        try:
            mutex_path = os.path.join(path, 'tcga-mutex-results')
        except Exception as e:
            raise BioagentException.PathNotFoundException()


        folders = os.listdir(mutex_path)
        with self.cadb:
            cur = self.cadb.cursor()
            cur.execute("DROP TABLE IF EXISTS Mutex")
            cur.execute("CREATE TABLE Mutex(Disease TEXT, Id1 TEXT, Id2 TEXT, Id3 TEXT, Id4, Id5, Score REAL)")
            for folder in folders:
                try:
                    if folder in tcga_study_names:
                        disease_path = os.path.join(mutex_path, folder)
                        file_path = os.path.join(disease_path, 'whole/no-network/ranked-groups.txt')
                    else:
                        continue
                except Exception as e:
                    raise BioagentException.PathNotFoundException()

                mutex_file = open(file_path, 'r')
                next(mutex_file)  # skip the header line
                for line in mutex_file:
                    vals = line.split('\t')
                    vals[len(vals) - 1] = vals[len(vals) - 1].rstrip('\n')
                    score = float(vals[0])
                    if score > 0.05:
                        continue

                    genes = []
                    for i in range(2, len(vals)):
                        genes.append(vals[i])
                        # print(genes[i-2])

                    # fill the rest with none
                    for i in range(len(vals), 7):
                        genes.append(None)


                    cur.execute("INSERT INTO Mutex VALUES(?, ?, ?, ?, ?,?, ?)",
                                (folder, genes[0], genes[1], genes[2], genes[3], genes[4], score))

        mutex_file.close()
    def populate_explained_table(self):
        """
        Find the correlations with a causal explanation
        :return:
        """
        with self.cadb:
            cur = self.cadb.cursor()
            cur.execute("DROP TABLE IF EXISTS Explained_Correlations")
            cur.execute("CREATE TABLE Explained_Correlations AS SELECT * FROM Correlations "
                        "LEFT JOIN Causality ON Causality.Id1 = Correlations.Id1 AND Causality.Id2 = Correlations.Id2  "
                        "AND Causality.PSite1 = Correlations.PSite1 AND Causality.PSite2 = Correlations.PSite2 "
                        "WHERE Rel IS NOT NULL",
                        ).fetchall()

    def populate_unexplained_table(self):
        """
        Find the correlations without a causal explanation
        :return:
        """
        with self.cadb:
            cur = self.cadb.cursor()
            cur.execute("DROP TABLE IF EXISTS Unexplained_Correlations")
            cur.execute("CREATE TABLE Unexplained_Correlations AS SELECT * FROM Correlations "
                        "LEFT JOIN Causality ON Causality.Id1 = Correlations.Id1 AND Causality.Id2 = Correlations.Id2  "
                        "AND Causality.PSite1 = Correlations.PSite1 AND Causality.PSite2 = Correlations.PSite2 "
                        "WHERE Rel IS NULL",
                        ).fetchall()


    def populate_sif_relations_table(self, path):
        """
        All sif relations from PathwayCommons
        :param path: Path to the folder that keeps PC.sif
        :return:
        """

        try:
            pc_path = os.path.join(path, 'PC.sif')
        except Exception as e:
            raise BioagentException.PathNotFoundException()

        pc_file = open(pc_path, 'r')

        with self.cadb:
            cur = self.cadb.cursor()
            cur.execute("DROP TABLE IF EXISTS Sif_Relations")
            cur.execute("CREATE TABLE Sif_Relations(Id1 TEXT,  Id2 TEXT, Rel TEXT)")
            for line in pc_file:
                vals = line.split('\t')
                id1 = vals[0].upper()
                id2 = (vals[2].rstrip('\n')).upper()
                rel = vals[1]
                cur.execute("INSERT INTO Sif_Relations VALUES(?, ?, ?)", (id1, id2, rel))

        pc_file.close()



    def populate_tcga_names_table(self, path):
        """
        Fills the mutation significance table for all genes and TCGA studies
        :param path: Path to the folder that keeps  tcga_disease_names.tsv
        :return:
        """

        try:
            tcga_path = os.path.join(path, 'tcga_disease_names.tsv')
        except Exception as e:
            raise BioagentException.PathNotFoundException()

        with open(tcga_path, 'rU') as tcga_file:
            tcga_file = csv.reader(tcga_file, delimiter='\t')

            with self.cadb:
                cur = self.cadb.cursor()
                cur.execute("DROP TABLE IF EXISTS TCGA")
                cur.execute("CREATE TABLE TCGA(LongName TEXT, Abbr TEXT)")

                next(tcga_file)  # skip the header line

                for row in tcga_file:
                    cur.execute("INSERT INTO TCGA VALUES(?, ?)",
                                (str(row[0]).lower(), str(row[1].rstrip('\n'))))


