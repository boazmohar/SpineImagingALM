# Encoding: utf-8
""" Step class
 Step is used to define an interface for transforming Thunder Imaging objects
"""
from prep.Steps import *


class Step(object):
    """ Step class defines an interface for transforming Thunder Imaging objects
        Has one method: do()
    """
    def __init__(self, func, **kwargs):
        self.func = func
        self.args = kwargs

    def __repr__(self):
        return 'Step: func: %s, kwargs: %s' % (self.func, self.args)

    def do(self, data):
        """ evaluates the function of ths step and passes on the given arguments

        :param data: Thunder Images object
        :return: Modified Thunder Images object
        """

        return eval(self.func)(data=data, **self.args)


