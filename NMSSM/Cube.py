#########################################################################
#                                                                       #
#    C u b e                                                            #
#                                                                       #
#########################################################################

# A class for saving an item in the cube, and keeping track of the variables in
# the cube by writing an information file.

# External modules.
import os
import sys

# Decorators for tracking cube etc.


class countcalls(object):

    """"Decorator that keeps track of the number of times a function is called.
    """
    __instances = {}

    def __init__(self, f):
        self.__f = f
        # Start count at -1, so that first call is zeroth.
        self.__numcalls = -1
        countcalls.__instances[f] = self

    def __call__(self, *args, **kwargs):
        "Function called - increment the counter."
        self.__numcalls += 1
        return self.__f(*args, **kwargs)

    def count(self):
        "Return the number of times the function f was called."
        return countcalls.__instances[self.__f].__numcalls

    def reset(self):
        "Reset the counter."
        # Start count at -1, so that first call is zeroth.
        self.__numcalls = -1


class memoize(object):

    """ Memoize a class's arguments and type. Though these
    are properties of the class, rather than the instance.
    """
    __instances = {}

    def __init__(self, f):
        self.__f = f
        self.__f.args = []
        self.__f.type = f

    def __call__(self, *args, **kwargs):
        self.__f.args = args
        return self.__f(*args, **kwargs)


class once(object):

    """ A function with this decorator is called only once.
    """
    __instances = {}

    def __init__(self, f):
        self.__f = f
        self.first = True

    def __call__(self, *args, **kwargs):
        if self.first:
            self.first = False
            return self.__f(*args, **kwargs)
        else:
            return


# List of labels for cube.
label = {}


@countcalls
def AddCube(cube, x, label, name):
    """ Function for adding the next item to the cube,
    and to its list of labels.
    Arguments:
    cube -- MultiNest cube to which value is appended.
    x -- value to add to cube.
    label -- Labels to which new label is appended.
    name -- New label to add to labels.
    """
    # Import variable inside function, to avoid circular module
    # imports.
    import MNOptions as MN

    # +1 because our count begins at zero.
    if AddCube.count() + 1 > MN.n_params:
        sys.exit(
            '''Cannot add more parameters to the cube than the length of the MultiNest cube.
                Have you added additonal parameters to the cube within telling MultiNest?
                i) Within MNOptions.py, increase n_params.''')

    cube[AddCube.count()] = x
    label[AddCube.count()] = name

@once
def PrintInfo():
    """ Function to write an info file, for use with plotting,
    and for saving Multinest options, likelihood functions and
    priors.
    """
    # Import variable inside function, to avoid circular module
    # imports.
    import MNOptions as MN
    import Priors
    import Likelihood

    infofile = MN.outputfiles_basename + '.info'
    info = open(infofile, 'w')

    info.write(
        '''# SuperPy info file, with cube indicies.
# This info file can be used for plot labels.\n''')

    for key, value in label.iteritems():
        # Add 1 to the keys, so that labelling begins at 1, rather than 0.
        # This is identical to the SuperBayeS info file convention.
        key = str(int(key) + 1)
        info.write('''lab%s=%s\n''' % (key, value))

    # Write the prior types and arguements.
    info.write('''# Priors.\n''')
    for key in Priors.CNMSSMModelTracker().param:
        info.write(
            '''%s %s %s\n''' %
            (key,
             Priors.CNMSSMModelTracker().param[key].type,
             Priors.CNMSSMModelTracker().param[key].arg))

    # Write the likelihoods types and arguments.
    info.write('''# Likelihoods.\n''')
    for key in Likelihood.CNMSSMConstraintTracker().constraint:
        info.write(
            '''%s %s %s\n''' %
            (key,
             Likelihood.CNMSSMConstraintTracker().constraint[key].type,
             Likelihood.CNMSSMConstraintTracker().constraint[key].arg))

    # Write the MN parameters.
    info.write('''# MultiNest options.\n''')
    for vars in dir(MN):
        if not vars.startswith("__"):
            info.write('''%s %s\n''' % (vars,  getattr(MN, vars)))
    info.close()
