import sys
import os
import json
import logging
from bioagents import Bioagent
from causality_agent import CausalityAgent
from indra.trips.processor import TripsProcessor

from kqml import KQMLModule, KQMLPerformative, KQMLList, KQMLString, KQMLToken
from bioagents.mra import MRA, MRA_Module
from bioagents.mra.mra_module import ekb_from_agent, get_target

logging.basicConfig(format='%(levelname)s: %(name)s - %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('CausalA')

_resource_dir = os.path.dirname(os.path.realpath(__file__)) + '/resources/'


class CausalityModule(KQMLModule):
    name = 'CausalA'
    def __init__(self, **kwargs):
        self.tasks = ['FIND-CAUSAL-PATH', 'FIND-CAUSALITY-TARGET', 'FIND-CAUSALITY-SOURCE',
                    'DATASET-CORRELATED-ENTITY', 'FIND-COMMON-UPSTREAMS',
                    'RESTART-CAUSALITY-INDICES']
        self.CA = CausalityAgent(_resource_dir)
        # Call the constructor of KQMLModule
        super(CausalityModule, self).__init__(**kwargs)
        # # Send subscribe messages
        # for task in self.tasks:
        #     self.subscribe_request(task)
        # # Send ready message
        # self.ready()
        # self.start()

    def receive_request(self, msg, content):
        """Handle request messages and respond.

        If a "request" message is received, decode the task and the content
        and call the appropriate function to prepare the response. A reply
        message is then sent back.
        """
        task_str = content.head().upper()
        if task_str == 'FIND-CAUSAL-PATH':
            try:
                reply_content = self.respond_find_causal_path(content)
            except Exception as e:
                logger.error(e)
                reply_content = make_failure()
        elif task_str == 'FIND-CAUSALITY_TARGET':
            try:
                reply_content = self.find_causality_target(content)
            except Exception as e:
                logger.error(e)
                reply_content = make_failure()
        elif task_str == 'FIND-CAUSALITY_SOURCE':
            try:
                reply_content = self.find_causality_source(content)
            except Exception as e:
                logger.error(e)
                reply_content = make_failure()
        elif task_str == 'DATASET-CORRELATED-ENTITY':
            try:
                reply_content = self.dataset_correlated_entity(content)
            except Exception as e:
                logger.error(e)
                reply_content = make_failure()
        elif task_str == 'FIND-COMMON-UPSTREAMS':
            try:
                reply_content = self.find_common_upstreams(content)
            except Exception as e:
                logger.error(e)
                reply_content = make_failure()
        elif task_str == 'RESTART-CAUSALITY-INDICES':
            try:
                reply_content = self.restart_causality_indices(content)
            except Exception as e:
                logger.error(e)
                reply_content = make_failure()
        else:
            reply_content = make_failure()

        reply_msg = KQMLPerformative('reply')
        reply_msg.set('content', reply_content)
        self.reply(msg, reply_msg)

    def make_indra_json(self, causality):
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
                indra_json = {'type': rel_type, 'enz': {'name': causality.id2,
                                                       'mods': [{'mod_type': 'phosphorylation', 'is_modified': True,
                                                                 'residue': causality.res2, 'position': causality.pos2}]},
                              'sub': {'name': causality.id1}, 'residue': causality.res1, 'position': causality.pos1}
            else:
                indra_json = {'type': rel_type, 'enz': {'name': causality.id1,
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

    def respond_find_causal_path(self, content):
        """Response content to find-qca-path request"""
        source_arg = content.get('SOURCE')
        target_arg = content.get('TARGET')

        if not source_arg:
            raise ValueError("Source is empty")
        if not target_arg:
            raise ValueError("Target is empty")

        print(source_arg)
        target = {'id': self._get_term_name(target_arg), 'pSite': ''}
        source = {'id': self._get_term_name(source_arg), 'pSite': ''}

        result = self.CA.find_causality({'source': source, 'target': target})


        if not result:
            reply = make_failure('NO_PATH_FOUND')
            return reply

        indra_json = [self.make_indra_json(result)]

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

def make_failure(reason=None):
    msg = KQMLList('FAILURE')
    if reason:
        msg.set('reason', reason)
    return msg

if __name__ == "__main__":
    CausalityModule(argv=sys.argv[1:])

cm = CausalityModule()
mapk1 = Agent('MAPK1', db_refs={'HGNC': '3236', 'TEXT': 'EGFR'})
term1 = ekb_from_agent(mapk1)
mapk2 = Agent('JUND', db_refs={'HGNC': '3236', 'TEXT': 'EGFR'})
term2 = ekb_from_agent(mapk1)
cm.respond_find_causal_path({'SOURCE': term1, 'TARGET': term2})