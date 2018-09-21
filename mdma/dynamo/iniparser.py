#!/usr/bin/env python3

import ast
try: from configparser import ConfigParser, RawConfigParser
except ImportError: from ConfigParser import ConfigParser, RawConfigParser

class ConfigFile(dict):
    def __init__(self, path=None):
        if path is not None: self.read(path)

    def read(self, path):
        config = ConfigParser()
        config.read(path)

        config = dict(config._sections)
        for section in config:
            if section not in self: self[section] = dict()
            for parameter in config[section]:
                value = ast.literal_eval(config[section][parameter])
                self[section][parameter] = value

    def write(self, path):
        config = RawConfigParser()

        for section in self:
            config.add_section(section)
            for parameter in self[section]:
                config.set(section, parameter, self[section][parameter])

        with open(path, 'w') as f: config.write(f)
