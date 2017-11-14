import re
from database_initializer import DatabaseInitializer


class CausalityAgent:
    def __init__(self, path):
        self.corr_ind = 0
        self.causality_ind = 0

        self.db_initializer = DatabaseInitializer(path)

        self.cadb = self.db_initializer.cadb

    def __del__(self):
        self.cadb.close()

    def reset_indices(self):
        self.corr_ind = 0
        self.causality_ind = 0

    def get_tcga_abbr(self, long_name):
        """
        Gets the study abbreviation given its long name
        :param long_name:
        :return:
        """
        with self.cadb:
            cur = self.cadb.cursor()

            name = cur.execute("SELECT Abbr FROM TCGA WHERE longName = ?", (long_name,)).fetchone()

            if not name:
                return None

            return name[0]


    @staticmethod
    def row_to_causality(row):
        """
          Convertd a row from sql table into causality object
          Positions need to be trimmed to correct PC formatting. E.g. s100S for pSite
        """

        sites1 = re.findall('([TYS][0-9]+)[TYS]', row[1])
        sites2 = re.findall('([TYS][0-9]+)[TYS]', row[3])
        mods1 = [{'mod_type': 'phosphorylation',
                  'residue': site[0],
                  'position': site[1:],
                  'is_modified': True}
                  for site in sites1]
        mods2 = [{'mod_type': 'phosphorylation',
                  'residue': site[0],
                  'position': site[1:],
                  'is_modified': True}
                  for site in sites2]
        causality = {'id1': row[0], 'mods1': mods1,
                     'id2': row[2], 'mods2': mods2,
                     'rel': row[4],
                     'uri_str': row[5]
                     }
        return causality


    @staticmethod
    def row_to_correlation(row):
        """
        Converts a row from sql table into correlation object
        Positions need to be trimmed to correct PC formatting. E.g. s100S for pSite
        :param row:
        :return:
        """
        l1 = len(row[1])
        p_site1 = row[1][1:l1-1]

        l2 = len(row[3])
        p_site2 = row[3][1:l2-1]

        corr = {'id1': row[0], 'pSite1': p_site1,
                'id2': row[2], 'pSite2': p_site2, 'correlation': row[4],
                'pVal': row[5],
                'explainable': "unassigned"}
        return corr

    def find_causality(self, param):
        """
        Finds the causal relationship between gene1 and gene2
        :param param: {source:{id: }, target:{id:}}
        :return:
        """

        with self.cadb:
            cur = self.cadb.cursor()
            sources = param.get('source').get('id')
            targets = param.get('target').get('id')

            if isinstance(sources, list):
                source_str = "(" + ", ".join((("'" + str(source) + "'") for source in sources)) + ")"
            else:
                source_str = "('" + sources + "')"

            if isinstance(targets, list):
                target_str = "(" + ", ".join((("'" + str(target) + "'") for target in targets)) + ")"
            else:
                target_str = "('" + targets + "')"

            query = "SELECT * FROM Causality WHERE Id1 IN " + source_str + "AND Id2 IN  " + target_str

            rows = cur.execute(query).fetchall()

            if len(rows) > 0:
                row = rows[0]
                causality = self.row_to_causality(row)
                return causality
            else:
                return ''

    def find_causality_targets(self, param):
        """
        Finds the causal relationship from gene list
        :param param: param: {id:[]}
        :return:
        """
        with self.cadb:
            cur = self.cadb.cursor()
            genes = param.get('id')

            if isinstance(genes, list):
                id_str =  ", ".join((("'" + str(gene) + "'") for gene in genes))
            else:
                id_str = "('" + genes + "')"

            rel = param.get('rel')

            if rel.upper() == "MODULATES":
                query = "SELECT * FROM Causality WHERE Id1 IN " + "(" + id_str + ")";
                rows = cur.execute(query).fetchall()
            else:
                query = "SELECT * FROM Causality WHERE Rel = ?  AND Id1 IN " + "(" + id_str + ")";
                rows = cur.execute(query, (rel,)).fetchall()

            targets = []
            for row in rows:
                causality = self.row_to_causality(row)
                targets.append(causality)

            return targets

    def find_next_correlation(self, gene):
        """
        Returns the next interesting relationship about gene. Can be explained or unexplained
        :param gene:
        :return:
        """

        with self.cadb:
            cur = self.cadb.cursor()

            causal_rows = cur.execute("SELECT * FROM Explained_Correlations WHERE Id1 = ? OR Id2 = ? "
                                      "ORDER BY ABS(Corr) DESC",
                                      (gene, gene)).fetchall()

            row_cnt = len(causal_rows)
            if row_cnt > self.causality_ind:
                row = causal_rows[self.causality_ind]
                self.causality_ind = self.causality_ind + 1

                corr = self.get_correlation_between(row[0], row[1], row[2], row[3])
                corr['explainable'] = "explainable"
            else:
                corr = self.find_next_unexplained_correlation(gene)

            # revert correlation info
            if corr != '' and corr['id2'] == gene:
                tmp = corr['id1']
                corr['id1'] = corr['id2']
                corr['id2'] = tmp
                tmp = corr['pSite1']
                corr['pSite1'] = corr['pSite2']
                corr['pSite2'] = tmp

            return corr

    def get_correlation_between(self, gene1, p_site1, gene2, p_site2):
        """
        When We are sure that there is a correlation between these
        :param gene1:
        :param p_site1:
        :param gene2:
        :param p_site2:
        :return:
        """
        with self.cadb:
            cur = self.cadb.cursor()
            # Don't change the order
            rows = cur.execute("SELECT * FROM Correlations WHERE Id1 = ? AND PSite1 = ?  AND Id2 = ?  AND PSite2 = ? "
                               "OR Id1 = ? AND PSite1 = ?  AND Id2 = ?  AND PSite2 = ? ",
                               (gene1, p_site1, gene2, p_site2, gene2, p_site2, gene1, p_site1)).fetchall()

            corr = ''
            if len(rows) > 0:
                row = rows[0]
                corr = self.row_to_correlation(row)

            return corr

    def find_next_unexplained_correlation(self, gene):
        """
        Finds the next highest unexplained correlation
        :param gene:
        :return:
        """
        with self.cadb:
            cur = self.cadb.cursor()
            rows = cur.execute("SELECT * FROM Unexplained_Correlations "
                               "WHERE Id1 = ? OR Id2 = ? ORDER BY ABS(Corr) DESC",
                               (gene, gene)).fetchall()

            row_cnt = len(rows)
            if row_cnt > self.corr_ind:
                row = rows[self.corr_ind]
                self.corr_ind = self.corr_ind + 1
                corr = self.row_to_correlation(row)
                corr['explainable'] = "unexplainable"
                return corr
            else:
                return ''

    def find_mutation_significance(self, gene, disease):
        """
        :param gene: single gene name
        :return: string, mutation significance
        """
        with self.cadb:
            cur = self.cadb.cursor()

            p_val = cur.execute("SELECT PVal FROM MutSig WHERE Id = ? AND Disease = ?", (gene, disease)).fetchone()

            if p_val[0] < 0.01:
                return 'highly significant'
            elif p_val[0] < 0.05:
                return 'significant'
            else:
                return 'not significant'

    def find_common_upstreams(self, genes):
        """
        Find common upstreams between a list of genes
        :param genes:
        :return:
        """

        with self.cadb:
            cur = self.cadb.cursor()

            if len(genes) < 2:
                return ''

            gene1 = genes[0]
            gene2 = genes[1]

            upstreams = cur.execute("SELECT s1.Id1 FROM Sif_Relations s1 "
                                    "INNER JOIN Sif_Relations s2 ON (s2.Id1 = s1.Id1 AND s1.Id2 = ? AND s2.id2 = ? AND  "
                                    "s1.Rel = 'controls-state-change-of' AND s2.Rel = s1.Rel)",
                                    (gene1, gene2)).fetchall()

            for i in range(2, len(genes)):
                gene = genes[i]

                upstream_arr = []
                for upstream in upstreams:
                    upstream_arr.append(upstream[0])

                query = "SELECT Id1 FROM Sif_Relations WHERE Rel = 'controls-state-change-of' AND Id2 = ? AND Id1 IN (" \
                        + ", ".join((("'"+str(upstream)+"'") for upstream in upstream_arr)) + ")"

                upstreams = cur.execute(query, (gene,)).fetchall()

            #format upstreams
            upstream_list = []
            for genes in upstreams:
                upstream_list.append(str(genes[0]).encode('utf8'))

            return upstream_list

    def find_mutex(self, gene):
        """Find a mutually exclusive group that includes gene
        ;:param single gene name
        :return object list
        """

        with self.cadb:
            cur = self.cadb.cursor()
            groups = cur.execute("SELECT * FROM Mutex WHERE Id1 = ? OR Id2 = ? OR Id3 = ?",
                                 (gene, gene, gene)).fetchall()

        # format groups
        mutex_list = []
        for group in groups:
            mutex = {'group': [], 'score': str(round(group[len(group) - 1], 2))}
            for i in range(len(group) - 1):
                if group[i] is not None:
                    mutex['group'].append(group[i].encode('utf8'))
            mutex_list.append(mutex)

        return mutex_list
#
# ca = CausalityAgent('./resources')
#
# ca.db_initializer.populate_tcga_names_table('./resources')

# print(ca.get_tcga_abbr("Ovarian serous cystadenocarcinoma"))
# res = ca.find_causality({'source': {'id':'MAPK1'}, 'target': {'id': ['JUND', 'ERF']}})

# uri_str = res['uri_str']

# pc_url = 'http://www.pathwaycommons.org/pc2/get?' + uri_str + 'format=SBGN'
# html = "<a href= \"" + pc_url + "\" target= \"_blank\">" + uri_str + "</a>"
# print(html)


# ca.find_causality_targets({'id': ['MAPK1', 'BRAF'], 'rel': 'phosphorylates'})

# ca.find_causality_targets({'id':'MAPK1', 'rel': 'phosphorylates'}, print_result)
# ca.find_causality_targets({'id':'BRAF', 'rel': 'is-phosphorylated-by'}, print_result)
#
# ca.populate_mutsig_table('./resources')
# ca.populate_sif_relations_table()

# ca.populate_explained_table()
# print(ca.find_next_unexplained_correlation('AKT1'))
# ca.find_next_unexplained_correlation('AKT1', print_result)
# ca.find_next_unexplained_correlation('AKT1', print_result)
# # ca.find_causality_targets({'source':{'id' :'BRAF', 'pSite' :'S365S', 'rel': 'is-phosphorylated-by'}},print_result)
# ca.find_next_correlation('AKT1',print_result)
# ca.find_next_correlation('AKT1',print_result)
# ca.find_next_correlation('AKT1',print_result)
# ca.find_next_correlation('AKT1',print_result)
# ca.find_next_correlation('AKT1',print_result)
# # ca.find_next_correlation('AKT1',print_result)
# ca.find_correlation_between('AKT1', 'BRAF')
# print(ca.find_mutation_significance('TP53', 'OV'))
# ca.find_common_upstreams('RAC1', 'RAC2')
# print(ca.find_common_upstreams(['AKT1', 'BRAF', 'MAPK1']))
