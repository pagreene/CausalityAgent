import re
import os
import logging
import sqlite3
import numpy as np
import socket
import logging
import uuid
from socketIO_client import SocketIO
import database_handler


class CausalityAgent:

    def __init__(self, sbgnviz_port = 3000):
        self.sbgnviz_port = sbgnviz_port
        self.user_id = '%s' % uuid.uuid4()
        self.user_name = 'CA'
        self.color_code = '#ff46a7'
        self.room_id = ''
        self.current_users = []

        self.db_handler = database_handler.DatabaseHandler()
        self.connect_sbgnviz()

    def __del__(self):
        self.db_handler.__del__

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

    def on_subscribe(self, room):
        event = 'subscribeAgent'
        self.room_id = room

        print('Reconnecting ' + room)

        user_info = {'userName': self.user_name,
                     'room': self.room_id,
                     'userId': self.user_id,
                     'colorCode': self.color_code}

        # self.socket_s.on('message', self.on_sbgnviz_message)
        self.socket_s.on('findCausality', self.db_handler.find_causality)
        self.socket_s.on('findCausalityTargets', self.db_handler.find_causality_targets)
        self.socket_s.on('findCorrelation', self.db_handler.find_next_correlation)
        self.socket_s.on('reconnect', self.connect_sbgnviz)
        self.socket_s.emit(event, user_info)
        self.socket_s.emit('agentNewFileRequest', {'room': self.room_id})
        self.socket_s.emit('agentConnectToTripsRequest', user_info)

    def on_user_list(self, user_list):
        self.current_users = user_list

    # def on_sbgnviz_message(self, data):
    #     if not isinstance(data, dict):
    #         return
    #     comment = data.get('comment')
    #     if isinstance(comment, list):
    #         comment = comment[0]
    #     if isinstance(comment, dict):
    #         comment = comment.get('text')
    #
    #     print(comment)


agent = CausalityAgent()
agent.start()

