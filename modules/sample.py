import os
import sys

class Sample():

    def __init__(self, name):
        self._name = name

    def add(self, key, value):
        '''
            Add a new attribute dinamically
        '''
        if not key or not value:
            raise ValueError("Expected key and val to be provided.")
        setattr(self, key, value)
    @property
    def name(self):
        return self._name
