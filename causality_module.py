import sys
import os
import json
import logging
from bioagents import Bioagent
from causality_agent import CausalityAgent
from indra.sources.trips.processor import TripsProcessor

from kqml import KQMLModule, KQMLPerformative, KQMLList, KQMLString, KQMLToken
from bioagents.mra import MRA, MRA_Module
from bioagents.mra.mra_module import ekb_from_agent, get_target

logging.basicConfig(format='%(levelname)s: %(name)s - %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('CausalA')

_resource_dir = os.path.dirname(os.path.realpath(__file__)) + '/resources/'


class CausalityModule(Bioagent):
    name = 'CausalA'
    tasks = ['FIND-CAUSAL-PATH', 'FIND-CAUSALITY-TARGET',
             'FIND-CAUSALITY-SOURCE',
             'DATASET-CORRELATED-ENTITY', 'FIND-COMMON-UPSTREAMS',
             'RESTART-CAUSALITY-INDICES']
    def __init__(self, **kwargs):
        self.CA = CausalityAgent(_resource_dir)
        # Call the constructor of KQMLModule
        super(CausalityModule, self).__init__(**kwargs)


    def respond_find_causal_path(self, content):
        """Response content to find-qca-path request"""
        source_arg = content.get('SOURCE')
        target_arg = content.get('TARGET')

        if not source_arg:
            raise ValueError("Source is empty")
        if not target_arg:
            raise ValueError("Target is empty")

        target = {'id': self._get_term_name(target_arg), 'pSite': ''}
        source = {'id': self._get_term_name(source_arg), 'pSite': ''}

        result = self.CA.find_causality({'source': source, 'target': target})


        if not result:
            reply = self.make_failure('NO_PATH_FOUND')
            return reply

        indra_json = [make_indra_json(result)]

        reply = KQMLList('SUCCESS')
        reply.set('paths', indra_json)

        return reply


    def _get_term_name(self, term_str):
        term_str = '<ekb>' + term_str + '</ekb>'
        tp = TripsProcessor(term_str)
        terms = tp.tree.findall('TERM')
        term_id = terms[0].attrib['id']
        agent = tp._get_agent_by_id(term_id, None)
        return agent.name

def make_indra_json(causality):
    """Convert causality response to indra format
        Causality format is (id1, res1, pos1, id2,res2, pos2, rel)"""
    if causality.pos1 == '':
        causality.pos1 = None
    if causality.pos2 == '':
        causality.pos2 = None
    if causality.res1 == '':
        causality.res1 = None
    if causality.res2 == '':
        causality.res2 = None

    causality.rel = causality.rel.upper()

    indra_relation_map = {
        "PHOSPHORYLATES": "Phosphorylation",
        "IS-PHOSPHORYLATED-BY": "Phosphorylation",
        "IS-DEPHOSPHORYLATED-BY": "Dephosphorylation",
        "UPREGULATES-EXPRESSION": "IncreaseAmount",
        "EXPRESSION-IS-UPREGULATED-BY": "IncreaseAmount",
        "DOWNREGULATES-EXPRESSION": "DecreaseAmount",
        "EXPRESSION-IS-DOWNREGULATED-BY": "DecreaseAmount"
    }

    rel_type = indra_relation_map[causality.rel]

    if "PHOSPHO" in causality.rel:  # phosphorylation
        if "IS" in causality.rel:  # passive
            indra_json = {'type': rel_type, 'enz':
                          {'name': causality.id2,
                            'mods': [{'mod_type': 'phosphorylation',
                                      'is_modified': True,
                                      'residue': causality.res2,
                                      'position': causality.pos2}]},
                          'sub': {'name': causality.id1},
                          'residue': causality.res1,
                          'position': causality.pos1}
        else:
            indra_json = {'type': rel_type,
                          'enz': {'name': causality.id1,
                                                   'mods': [{'mod_type': 'phosphorylation', 'is_modified': True,
                                                             'residue': causality.res1,
                                                             'position': causality.pos1}]},
                          'sub': {'name': causality.id2}, 'residue': causality.res2, 'position': causality.pos2}
    else:  # regulation
        if "IS" in causality.rel:  # passive
            indra_json = {'type': rel_type, 'subj': {'name': causality.id2,
                                                   'mods': [{'mod_type': 'phosphorylation', 'is_modified': True,
                                                             'residue': causality.res2, 'position': causality.pos2}]},
                          'obj': {'name': causality.id1}, 'residue': causality.res1, 'position': causality.pos1}
        else:
            indra_json = {'type': rel_type, 'subj': {'name': causality.id1,
                                                   'mods': [{'mod_type': 'phosphorylation', 'is_modified': True,
                                                             'residue': causality.res1,
                                                             'position': causality.pos1}]},
                          'obj': {'name': causality.id2}, 'residue': causality.res2, 'position': causality.pos2}

    return indra_json

if __name__ == "__main__":
    CausalityModule(argv=sys.argv[1:])

