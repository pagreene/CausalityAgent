import sys
import os
import json
import logging
from bioagents import Bioagent
from causality_agent import CausalityAgent
from indra.sources.trips.processor import TripsProcessor
from kqml import KQMLModule, KQMLPerformative, KQMLList, KQMLString, KQMLToken


logging.basicConfig(format='%(levelname)s: %(name)s - %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('CausalA')

_resource_dir = os.path.dirname(os.path.realpath(__file__)) + '/resources/'

cancerTypes = [
        {'abbr': "ACC", 'longName':"Adrenocortical carcinoma"},
        {'abbr': "BLCA",'longName':"Bladder Urothelial Carcinoma"},
        {'abbr':"BRCA",'longName':"Breast invasive carcinoma"},
        {'abbr':"CESC",'longName':"Cervical squamous cell carcinoma and endocervical adenocarcinoma"},
        {'abbr':"CHOL",'longName':"Cholangiocarcinoma"},
        {'abbr':"COAD",'longName':"Colon adenocarcinoma"},
        {'abbr':"COADREAD",'longName':"Colorectal cancer"},
        {'abbr':"DLBC",'longName':"Lymphoid Neoplasm Diffuse Large B-cell Lymphoma"},
        {'abbr':"GBM",'longName':"Glioblastoma multiforme"},
        {'abbr':"GBMLGG",'longName':"Glioblastoma multiforme/Brain Lower Grade Glioma"},
        {'abbr':"HNSC",'longName':"Head and Neck squamous cell carcinoma"},
        {'abbr':"KICH",'longName':"Kidney Chromophobe"},
        {'abbr':"KIPAN",'longName':"Pan-Kidney"},
        {'abbr':"KIRC",'longName':"Kidney renal clear cell carcinoma"},
        {'abbr':"KIRP",'longName':"Kidney renal papillary cell carcinoma"},
        {'abbr':"LAML",'longName':"	Acute Myeloid Leukemia"},
        {'abbr':"LGG",'longName':"Brain Lower Grade Glioma"},
        {'abbr':"LIHC",'longName':"Liver hepatocellular carcinoma"},
        {'abbr':"LUAD",'longName':"Lung adenocarcinoma"},
        {'abbr':"LUSC",'longName':"Lung squamous cell carcinoma"},
        {'abbr':"OV",'longName':"Ovarian serous cystadenocarcinoma"},
        {'abbr':"PAAD",'longName':"Pancreatic adenocarcinoma"},
        {'abbr':"PCPG",'longName':"Pheochromocytoma and Paraganglioma"},
        {'abbr':"PRAD",'longName':"Prostate adenocarcinoma"},
        {'abbr':"READ",'longName':"Rectum adenocarcinoma"},
        {'abbr':"SARC", 'longName':"Sarcoma"},
        {'abbr':"SKCM", 'longName':"Skin Cutaneous Melanoma"},
        {'abbr':"STAD", 'longName':"Stomach adenocarcinoma"},
        {'abbr':"STES", 'longName':"Stomach and Esophageal carcinoma"},
        {'abbr':"TGCT", 'longName':"Testicular Germ Cell Tumors"},
        {'abbr':"THCA", 'longName':"Thyroid carcinoma"},
        {'abbr':"UCEC", 'longName':"Uterine Corpus Endometrial Carcinoma"},
        {'abbr':"UCS", 'longName':"Uterine Carcinosarcoma"},
        {'abbr':"UVM", 'longName':"Uveal Melanoma"}]


class CausalityModule(Bioagent):
    name = 'CausalA'
    tasks = ['FIND-CAUSAL-PATH', 'FIND-CAUSALITY-TARGET',
             'FIND-CAUSALITY-SOURCE',
             'DATASET-CORRELATED-ENTITY', 'FIND-COMMON-UPSTREAMS',
             'RESTART-CAUSALITY-INDICES', 'FIND-MUTATION-SIGNIFICANCE', 'FIND-MUTEX']

    def __init__(self, **kwargs):
        self.CA = CausalityAgent(_resource_dir)
        # Call the constructor of KQMLModule
        super(CausalityModule, self).__init__(**kwargs)

    def respond_find_causal_path(self, content):
        """Response content to find-causal-path request"""

        source_arg = content.gets('SOURCE')
        target_arg = content.gets('TARGET')

        if not source_arg:
            raise ValueError("Source is empty")
        if not target_arg:
            raise ValueError("Target is empty")

        target_names = _get_term_names(target_arg)
        source_names = _get_term_names(source_arg)

        if not target_names or not source_names:
            reply = self.make_failure('NO_PATH_FOUND')
            return reply

        source_name = source_names[0]
        target_name = target_names[0]

        target = {'id': target_name, 'pSite': ''}
        source = {'id': source_name, 'pSite': ''}

        result = self.CA.find_causality({'source': source, 'target': target})

        if not result:
            reply = self.make_failure('NO_PATH_FOUND')
            return reply

        indra_json = json.dumps([make_indra_json(result)])

        reply = KQMLList('SUCCESS')
        reply.sets('paths', indra_json)

        # Send PC links to provenance tab
        self.send_provenance(result['uri_str'])

        return reply

    def send_provenance(self, uri_str):
        pc_url = 'http://www.pathwaycommons.org/pc2/get?' + uri_str + 'format=SBGN'
        html = '<a href= \'' + pc_url + '\' target= \'_blank\' > PC link</a>'
        msg = KQMLPerformative('tell')
        content = KQMLList('add-provenance')
        content.sets('html', html)
        pc_url_formatted = "http://www.pathwaycommons.org/pc2/get?" + uri_str + "format=SBGN"
        content.sets('pc', pc_url_formatted)
        msg.set('content', content)
        self.send(msg)

    def respond_find_causality_target(self, content):
        """Response content to find-causality-target request"""
        target_arg = content.gets('TARGET')
        rel = content.gets('TYPE')

        if not target_arg:
            raise ValueError("Target is empty")
        
        target_names = _get_term_names(target_arg)
        if not target_names:
            reply = self.make_failure('MISSING_MECHANISM')
            return reply
        target_name = target_names[0]


        rel_map = {
            "phosphorylation": "phosphorylates",
            "dephosphorylation": "dephosphorylates",
            "activate": "upregulates-expression",
            "increase": "upregulates-expression",
            "inhibit": "downregulates-expression",
            "decrease": "downregulates-expression",
            "modulate": "modulates",
        }

        target = {'id': target_name, 'pSite': ' ', 'rel': rel_map[rel]}
        result = self.CA.find_causality_targets(target)

        if not result:
            reply = self.make_failure('MISSING_MECHANISM')
            return reply

        # Send PC links to provenance tab
        # Multiple interactions are sent separately
        for r in result:
            self.send_provenance(r['uri_str'])

        indra_json = json.dumps([make_indra_json(r) for r in result])

        reply = KQMLList('SUCCESS')
        reply.sets('paths', indra_json)

        return reply

    def respond_find_causality_source(self, content):
        """Response content to find-qca-path request"""
        source_arg = content.gets('SOURCE')
        rel = content.gets('TYPE')

        if not source_arg:
            raise ValueError("Source is empty")

        source_names = _get_term_names(source_arg)
        if not source_names:
            reply = self.make_failure('MISSING_MECHANISM')
            return reply
        source_name = source_names[0]

        rel_map = {
            "phosphorylation": "is-phosphorylated-by",
            "dephosphorylation": "is-dephosphorylated-by",
            "activate": "expression-is-up",
            "increase": "upregulates-expression",
            "inhibit": "downregulates-expression",
            "decrease": "downregulates-expression",
            "modulate": "modulates",
        }

        source = {'id': source_name, 'pSite': ' ','rel': rel_map[rel]}

        result = self.CA.find_causality_targets(source)

        if not result:
            reply = self.make_failure('MISSING_MECHANISM')
            return reply

        # Multiple interactions are sent separately
        for r in result:
            self.send_provenance(r['uri_str'])

        indra_json = json.dumps([make_indra_json(r) for r in result])

        reply = KQMLList('SUCCESS')
        reply.sets('paths', indra_json)

        return reply

    def respond_dataset_correlated_entity(self, content):
        source_arg = content.gets('SOURCE')
        if not source_arg:
            raise ValueError("Source is empty")

        source_names = _get_term_names(source_arg)
        if not source_names:
            reply = self.make_failure('MISSING_MECHANISM')
            return reply
        source_name = source_names[0]
        res = self.CA.find_next_correlation(source_name)
        reply = KQMLList('SUCCESS')
        reply.sets('target', res['id2'])
        reply.sets('correlation', str(res['correlation']))
        reply.sets('explainable', res['explainable'])

        return reply;


    def respond_find_common_upstreams(self, content):
        """Response content to find-common-upstreams request"""

        genes_arg = content.gets('GENES')

        if not genes_arg:
            raise ValueError("Gene list is empty")

        gene_names = _get_term_names(genes_arg)

        if not gene_names:
            reply = self.make_failure('MISSING_MECHANISM')
            return reply

        gene_list = []
        for gene_name in gene_names:
            gene_list.append(str(gene_name))

        print(gene_list)
        result = self.CA.find_common_upstreams(gene_list)

        if not result:
            reply = self.make_failure('MISSING_MECHANISM')
            return reply

        reply = KQMLList('SUCCESS')
        reply.sets('upstreams', result)

        return reply

    def respond_find_mutation_significance(self, content):
        """Response content to find-mutation-significance request"""
        gene_arg = content.gets('GENE')

        disease_arg = content.gets('DISEASE')

        print("disease_arg")
        print(disease_arg)

        if not gene_arg:
            raise ValueError("Source is empty")

        gene_names = _get_term_names(gene_arg)
        if not gene_names:
            reply = self.make_failure('MISSING_MECHANISM')
            return reply
        gene_name = gene_names[0]

        if not disease_arg:
            result = self.CA.find_mutation_significance(gene_name, 'OV') # default is ovarian cancer
        else:
            disease_names = _get_term_names(disease_arg)
            if not disease_names:
                reply = self.make_failure('MISSING_MECHANISM')
                return reply
            disease_name = disease_names[0]
            print(disease_name)
            result = self.CA.find_mutation_significance(gene_name, disease_name)  # default is ovarian cancer

        if not result:
            reply = self.make_failure('MISSING_MECHANISM')
            return reply

        reply = KQMLList('SUCCESS')
        reply.sets('mutsig', result)

        return reply


    def respond_find_mutex(self, content):
        """Response content to find-mutation-significance request"""

        gene_arg = content.gets('GENE')

        if not gene_arg:
            raise ValueError("Source is empty")

        gene_names = _get_term_names(gene_arg)
        if not gene_names:
            reply = self.make_failure('MISSING_MECHANISM')
            return reply
        gene_name = gene_names[0]

        result = self.CA.find_mutex(gene_name)

        if not result:
            reply = self.make_failure('MISSING_MECHANISM')
            return reply

        reply = KQMLList('SUCCESS')
        reply.sets('mutex', result)

        return reply


def _get_term_names(term_str):
    """Given an ekb-xml returns the names of genes in a list"""

    tp = TripsProcessor(term_str)
    terms = tp.tree.findall('TERM')
    if not terms:
        return None

    agent_names = []
    for term in terms:
        term_id = term.attrib['id']
        agent = tp._get_agent_by_id(term_id, None)

        if agent is not None:
            if isinstance(agent, list):
                for a in agent:
                    if a.name:
                        agent_names.append(a.name)
            else:
                agent_names.append(agent.name)

    if len(agent_names) == 0:
        return None

    return agent_names

def make_indra_json(causality):
    """Convert causality response to indra format
        Causality format is (id1, res1, pos1, id2,res2, pos2, rel)"""

    causality['rel'] = causality['rel'].upper()

    indra_relation_map = {
        "PHOSPHORYLATES": "Phosphorylation",
        "IS-PHOSPHORYLATED-BY": "Phosphorylation",
        "IS-DEPHOSPHORYLATED-BY": "Dephosphorylation",
        "UPREGULATES-EXPRESSION": "IncreaseAmount",
        "EXPRESSION-IS-UPREGULATED-BY": "IncreaseAmount",
        "DOWNREGULATES-EXPRESSION": "DecreaseAmount",
        "EXPRESSION-IS-DOWNREGULATED-BY": "DecreaseAmount"
    }

    rel_type = indra_relation_map[causality['rel']]

    s, t = ('2', '1') if 'IS' in causality['rel'] else ('1', '2')
    subj, obj = ('enz', 'sub') if 'PHOSPHO' in causality['rel'] else \
                ('subj', 'obj')

    if "PHOSPHO" in causality['rel']:  # phosphorylation
        indra_json = {'type': rel_type,
                      subj: {'name': causality['id%s' % s],
                             'mods': causality['mods%s' % s]},
                      obj: {'name': causality['id%s' % t]},
                      'residue': causality['mods%s' % t][0]['residue'],
                      'position': causality['mods%s' % t][0]['position']}
    return indra_json




def make_indra_json(causality):
    """Convert causality response to indra format
        Causality format is (id1, res1, pos1, id2,res2, pos2, rel)"""

    causality['rel'] = causality['rel'].upper()

    indra_relation_map = {
        "PHOSPHORYLATES": "Phosphorylation",
        "IS-PHOSPHORYLATED-BY": "Phosphorylation",
        "IS-DEPHOSPHORYLATED-BY": "Dephosphorylation",
        "UPREGULATES-EXPRESSION": "IncreaseAmount",
        "EXPRESSION-IS-UPREGULATED-BY": "IncreaseAmount",
        "DOWNREGULATES-EXPRESSION": "DecreaseAmount",
        "EXPRESSION-IS-DOWNREGULATED-BY": "DecreaseAmount"
    }

    rel_type = indra_relation_map[causality['rel']]

    s, t = ('2', '1') if 'IS' in causality['rel'] else ('1', '2')
    subj, obj = ('enz', 'sub') if 'PHOSPHO' in causality['rel'] else \
                ('subj', 'obj')

    if "PHOSPHO" in causality['rel']:  # phosphorylation
        indra_json = {'type': rel_type,
                      subj: {'name': causality['id%s' % s],
                             'mods': causality['mods%s' % s]},
                      obj: {'name': causality['id%s' % t]},
                      'residue': causality['mods%s' % t][0]['residue'],
                      'position': causality['mods%s' % t][0]['position']}
    return indra_json


if __name__ == "__main__":
    CausalityModule(argv=sys.argv[1:])

