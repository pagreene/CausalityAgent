import logging
import sys
import uuid
from socketIO_client import SocketIO
import causality_agent
import os
import threading


_resource_dir = os.path.dirname(os.path.realpath(__file__)) + '/resources/'

class CausalitySbgnvizInterface(object):

    def __init__(self, sbgnviz_port = 3000):
        self.sbgnviz_port = sbgnviz_port
        self.user_id = '%s' % uuid.uuid4()
        self.user_name = 'CA'
        self.color_code = '#ff46a7'
        self.room_id = ''
        self.current_users = []
        if len(sys.argv) == 1:
            path = _resource_dir
        else:
            path = sys.argv[1]

        self.CA = causality_agent.CausalityAgent(path)
        self.connect_sbgnviz()


    def __del__(self):
        self.CA.__del__

    def start(self):
        # Wait for things to happen
        while True:
            try:
                self.socket_s.wait(seconds=0.1)
            except KeyboardInterrupt:
                break
        self.socket_s.emit('disconnect')
        self.socket_s.disconnect()

    def connect_sbgnviz(self):
        # Initialize sockets
        self.socket_s = SocketIO('localhost', self.sbgnviz_port)
        self.socket_s.emit('agentCurrentRoomRequest', self.on_subscribe)

    def send_connection_request(self):
        self.socket_s.emit('agentCurrentRoomRequest', self.on_subscribe)

    def on_subscribe(self, room):
        event = 'subscribeAgent'

        if room is None:
            threading.Timer(0.125, self.connect_sbgnviz).start()
            # return
        else:

            self.room_id = room

            user_info = {'userName': self.user_name,
                         'room': self.room_id,
                         'userId': self.user_id,
                         'colorCode': self.color_code}

            # self.socket_s.on('message', self.on_sbgnviz_message)
            self.socket_s.on('findCausality', self.on_find_causality)
            self.socket_s.on('findCausalityTargets',self.on_find_causality_targets)
            self.socket_s.on('findCorrelation', self.on_find_next_correlation)
            self.socket_s.on('findCommonUpstreams', self.on_find_common_upstreams)
            self.socket_s.on('reconnect', self.connect_sbgnviz)
            self.socket_s.emit(event, user_info)
            self.socket_s.emit('agentNewFileRequest', {'room': self.room_id})
            self.socket_s.emit('agentConnectToTripsRequest', user_info)

            print('Connected ' + room)

    def on_user_list(self, user_list):
        self.current_users = user_list

    def on_find_causality_targets(self, params, callback):
        res = self.CA.find_causality_targets(params)
        callback(res)

    def on_find_causality(self, params, callback):
        res = self.CA.find_causality
        callback(res)

    def on_find_next_correlation(self, params, callback):
        res = self.CA.find_next_correlation
        callback(res)

    def on_find_common_upstreams(self, params, callback):
        res = self.CA.find_common_upstreams
        callback(res)


if __name__ == '__main__':
    agent_interface = CausalitySbgnvizInterface()
    agent_interface.start()

