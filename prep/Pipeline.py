# Encoding: utf-8
""" Pipeline class
 Pipeline is a way to preform steps in order on a thunder Images object
"""
from __future__ import print_function
import sys
import collections
from prep.Step import Step


class Pipeline(object):
    """ Pipeline class aggregates Step classes in sequence

    """
    def __init__(self, source):
        self.steps = collections.OrderedDict()
        self.source = source

    def __repr__(self):
        if hasattr(self, 'source') and  hasattr(self, 'steps'):
            return 'Pipeline: from %s, steps: %s' % (self.source, self.steps)
        else:
            return 'Pipeline'

    def addStep(self, name, func, **kwargs):
        """ Adding a step to the ordered dict

        :param name: name of current step
        :param func: function name to call
        :param kwargs: arguments to give the function
        """
        self.steps[name] = Step(func, **kwargs)

    def doStep(self, data, name):
        """ evoke the do method of the step 'name' on data
        :param data: Images object
        :param name: name of step
        :return: Modified Images object
        """
        if name in self.steps:
            print('step:' + name)
            sys.stdout.flush()
            return self.steps[name].do(data)
        else:
            print('No step called: ' + name)

    def doAllSteps(self, data):
        """ evoke all step in pipeline in order
        :param data: Images object
        :return: Modified Images object
        """
        for k, v in self.steps.items():
            data = self.doStep(data, k)
        return data

    def remStep(self, name):
        """ remove a step from the dictionary
        :param name: the name of the step to remove
        """
        if name in self.steps:
            del self.steps[name]
        else:
            print('No step called: ' + name)
