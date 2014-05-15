#########################################################################
#                                                                       #
#    D e b u g                                                          #
#    Set debug options here.                                            #
#                                                                       #
#########################################################################

# External modules.
import faulthandler
import resource
import sys

#########################################################################

# When not in debug mode, this class in used to disable printing to the screen.


class NullWriter(object):

    ''' Disables print statements - useful for disabling lots of output. '''

    def write(self, arg):
        pass

#########################################################################

# If you want to diable printing to the screen and the seg-fault watcher,
# turn to False. This can be desirable, because the code can produce big log
# files.
Debug = True

if Debug:
    # Turn on seg fault handler.
    faulthandler.enable()

    # There is a note in the MultiNest documentation about increasing the stack
    # size to 256mb to avoid stack overflows or shifts that result in seg faults.
    # 262144 kb = 256 mb.
    resource.setrlimit(resource.RLIMIT_STACK, (262144, 262144))

else:
    # Disable printing to the screen from python.
    sys.stdout = NullWriter()
