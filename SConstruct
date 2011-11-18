# -*- python -*-
from lsst.sconsUtils import scripts, env

scripts.BasicSConstruct.initialize("ap")

# If one of the randomized unit tests fail, uncomment and
# set the following defines to obtain repeatable behaviour
#env.Append(CPPFLAGS = ' -DLSST_AP_INIT_RANDOM_SEC=')
#env.Append(CPPFLAGS = ' -DLSST_AP_INIT_RANDOM_NSEC=')

scripts.BasicSConstruct.finish()
