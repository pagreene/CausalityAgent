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
        source_arg = content.gets('SOURCE')
        target_arg = content.gets('TARGET')

        if not source_arg:
            raise ValueError("Source is empty")
        if not target_arg:
            raise ValueError("Target is empty")

        target_name = _get_term_name(target_arg)
        source_name = _get_term_name(source_arg)
        if not target_name or not source_name:
            reply = make_failure('NO_PATH_FOUND')
            return reply

        target = {'id': target_name, 'pSite': ''}
        source = {'id': source_name, 'pSite': ''}

        result = self.CA.find_causality({'source': source, 'target': target})


        if not result:
            reply = self.make_failure('NO_PATH_FOUND')
            return reply

        indra_json = json.dumps([make_indra_json(result)])

        reply = KQMLList('SUCCESS')
        reply.sets('paths', indra_json)

        return reply

    def respond_find_causality_target(self, content):
        """Response content to find-qca-path request"""
        target_arg = content.gets('TARGET')
        rel = content.gets('TYPE')

        if not target_arg:
            raise ValueError("Target is empty")

        target_name = _get_term_name(target_arg)
        if not target_name:
            reply = make_failure('MISSING_MECHANISM')
            return reply

        rel_map = {
            "phosphorylation": "phosphorylates",
            "dephosphorylation": "dephosphorylates",
            "activate": "upregulates-expression",
            "increase": "upregulates-expression",
            "inhibit": "downregulates-expression",
            "decrease": "downregulates-expression",
            "modulate": "modulates",
        }

        target = {'id': target_name, 'pSite': ' ',
                  'rel': rel_map[rel]}

        result = self.CA.find_causality_targets(target)


        if not result:
            reply = self.make_failure('MISSING_MECHANISM')
            return reply

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

        source_name = _get_term_name(source_arg)
        if not source_name:
            reply = make_failure('MISSING_MECHANISM')
            return reply

        rel_map = {
            "phosphorylation": "is-phosphorylated-by",
            "dephosphorylation": "is-dephosphorylated-by",
            "activate": "expression-is-upregulated-by",
            "increase": "expression-is-upregulated-by",
            "inhibit": "expression-is-downregulated-by",
            "decrease": "expression-is-downregulated-by",
            "modulate": "modulates",
        }


        source = {'id': source_name, 'pSite': ' ',
                  'rel': rel_map[rel]}

        result = self.CA.find_causality_targets(source)


        if not result:
            reply = self.make_failure('MISSING_MECHANISM')
            return reply

        indra_json = json.dumps([make_indra_json(r) for r in result])

        reply = KQMLList('SUCCESS')
        reply.sets('paths', indra_json)

        return reply

def _get_term_name(term_str):
    tp = TripsProcessor(term_str)
    terms = tp.tree.findall('TERM')
    if not terms:
        return None
    term_id = terms[0].attrib['id']
    agent = tp._get_agent_by_id(term_id, None)
    if agent is None:
        return None
    return agent.name

def make_indra_json(causality):
    """Convert causality response to indra format
        Causality format is (id1, res1, pos1, id2,res2, pos2, rel)"""
    # TODO: Do these special cases still need to be handled?

    '''
    if causality.pos1 == '':
        causality.pos1 = None
    if causality.pos2 == '':
        causality.pos2 = None
    if causality.res1 == '':
        causality.res1 = None
    if causality.res2 == '':
        causality.res2 = None
    '''

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


def make_failure(reason=None):
     msg = KQMLList('FAILURE')
     if reason:
         msg.set('reason', reason)
     return msg

if __name__ == "__main__":
    CausalityModule(argv=sys.argv[1:])

