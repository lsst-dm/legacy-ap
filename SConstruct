# -*- python -*-
import glob, os.path, re, os
import lsst.SConsUtils as scons

#
# Custom configure tests
#
visCheckSrc = """
    __attribute__((visibility('hidden')))  void hiddenFunc() {}
    __attribute__((visibility('default'))) void defaultFunc() {}
    int main(int argc, char **argv) {
        hiddenFunc();
        defaultFunc();
        return 0;
    }
    """

popcountCheckSrc = """
    int main(int argc, char **argv) {
        unsigned long long ull = 0;
        unsigned long      ul  = 0;
        unsigned int       ui  = 0;
        return __builtin_popcount(ui) + __builtin_popcountl(ul) + __builtin_popcountll(ull);
    }
    """

noatimeCheckSrc = """
    #include <sys/types.h>
    #include <sys/stat.h>
    #include <fcntl.h>
    int main(int argc, char **argv) {
        open('/tmp/dummy', O_RDONLY | O_NOATIME, 0);
        return 0;
    }
    """

nocacheCheckSrc = """
    #include <sys/types.h>
    #include <sys/stat.h>
    #include <fcntl.h>
    int main(int argc, char **argv) {
        fcntl(-1, F_NOCACHE, 1);
        return 0;
    }
    """

rshiftCheckSrc = """
    int main(int argc, char **argv) {
        char test[-1 >> 1];
        return 0;
    }
    """

long64CheckSrc = """
    int main(int argc, char **argv) {
        char test[sizeof(long) - 8];
        return 0;
    }
    """


def CustomCompilerFlag(context, flag):
    context.Message('Checking if compiler supports ' + flag + ' flag ')
    ccflagsOld = context.env['CCFLAGS']; 
    context.env.Append(CCFLAGS = flag)
    result = context.TryCompile("""
        int main(int argc, char **argv) { return 0; }
        """, '.c')
    context.Result(result)
    if not result:
        context.env.Replace(CCFLAGS = ccflagsOld)
    return result

def CustomCompileCheck(context, message, source, extension = '.c'):
    context.Message(message)
    result = context.TryCompile(source, extension)
    context.Result(result)
    return result

def CustomLinkCheck(context, message, source, extension = '.c'):
    context.Message(message)
    result = context.TryLink(source, extension)
    context.Result(result)
    return result

#
# Setup our environment
#
env = scons.makeEnv('associate',
                    r"$HeadURL$",
                    [['boost', 'boost/version.hpp',
                        'boost_filesystem boost_regex boost_serialization boost_program_options:C++'],
                     ["mysqlclient", "mysql/mysql.h", "mysqlclient:C++"],
                     ['python', 'Python.h'],
                     ['daf_base', 'lsst/daf/base/Citizen.h', 'daf_base:C++'],
                     ['pex_exceptions', 'lsst/pex/exceptions/Exception.h', 'pex_exceptions:C++'],
                     ['daf_data', 'lsst/daf/data/SupportFactory.h', 'daf_data:C++'],
                     ['afw', 'lsst/afw/detection/Source.h', 'afw:C++']
                    ])

#
# Run configure tests
#
if not env.CleanFlagIsSet():
    conf = Configure(env, custom_tests = {'CustomCompilerFlag' : CustomCompilerFlag,
                                          'CustomCompileCheck' : CustomCompileCheck,
                                          'CustomLinkCheck'    : CustomLinkCheck})
    if env['PLATFORM'] == 'posix':
        # POSIX platforms have AIO functionality in librt
        if not conf.CheckLibWithHeader('rt', 'aio.h', 'C'):
            print 'Missing support for Posix AIO'
            Exit(1)
        # Necessary to get Linux direct I/O support
        if os.uname()[0].title() == 'Linux':
            conf.env.Append(CPPFLAGS = ' -D_GNU_SOURCE')

    # If one of the randomized unit tests fail, uncomment and
    # set the following defines to obtain repeatable behaviour
    #conf.env.Append(CPPFLAGS = ' -DLSST_AP_INIT_RANDOM_SEC=')
    #conf.env.Append(CPPFLAGS = ' -DLSST_AP_INIT_RANDOM_NSEC=')

    # Required for [U]INT64_C()
    conf.env.Append(CPPFLAGS = ' -D__STDC_CONSTANT_MACROS')
    # Indicate that shared libraries are being built
    conf.env.Append(CPPFLAGS = ' -DLSST_AP_SHARED_LIBRARY_BUILD=1')
    # compiler flags and features
    conf.CustomCompilerFlag('-fvisibility-inlines-hidden')
    if conf.CustomCompileCheck('Checking for __attribute__((visibility)) support... ', visCheckSrc):
        conf.env.Append(CPPFLAGS = ' -DLSST_AP_HAVE_VISIBILITY=1')
    if conf.CustomCompileCheck('Checking for __builtin_popcount... ', popcountCheckSrc):
        conf.env.Append(CPPFLAGS = ' -DLSST_AP_HAVE_BUILTIN_POPCOUNT=1')
    if not conf.CustomCompileCheck('Checking for unsigned right shift ... ', rshiftCheckSrc):
        conf.env.Append(CPPFLAGS = ' -DLSST_AP_HAVE_SIGNED_RSHIFT=1')
    # Without some help, SWIG disagrees with boost on the actual type of int64_t
    if conf.CustomCompileCheck('Checking whether long is at least 8 bytes ... ', long64CheckSrc):
        conf.env.Append(SWIGFLAGS = '-DSWIGWORDSIZE64')
    # Platform features
    if conf.CheckFunc('clock_gettime'): # Linux/Solaris: prototype in <time.h>
        conf.env.Append(CPPFLAGS = ' -DLSST_AP_HAVE_CLOCK_GETTIME=1')
    if conf.CheckFunc('directio'): # Solaris: prototype in <sys/fcntl.h>
        conf.env.Append(CPPFLAGS = ' -DLSST_AP_HAVE_DIRECTIO=1')
    if conf.CustomCompileCheck('Checking for O_NOATIME in open()... ', noatimeCheckSrc):
        conf.env.Append(CPPFLAGS = ' -DLSST_AP_HAVE_O_NOATIME=1')
    if conf.CustomCompileCheck('Checking for F_NOCACHE in fcntl()... ', nocacheCheckSrc):
        conf.env.Append(CPPFLAGS = ' -DLSST_AP_HAVE_F_NOCACHE=1')
    env = conf.Finish()

#
# Build/install things
#
for d in Split('include/lsst/ap doc lib python/lsst/ap bin tests'):
    SConscript(os.path.join(d, 'SConscript'))

env['IgnoreFiles'] = r"(~$|\.pyc$|^\.svn$|\.o$)"

Alias('install', [env.Install(env['prefix'], 'python'),
                  env.Install(env['prefix'], 'include'),
                  env.Install(env['prefix'], 'bin'),
                  env.Install(env['prefix'], 'sql'),
                  env.Install(env['prefix'], 'lib'),
                  env.Install(env['prefix'] + '/doc', 'doc/html'),
                  env.InstallEups(env['prefix'] + '/ups', glob.glob('ups/*.table'))])

scons.CleanTree(r"*~ core *.so *.os *.o")

#
# Build TAGS files
#
files = scons.filesToTag()
if files:
    env.Command('TAGS', files, 'etags -o $TARGET $SOURCES')

env.Declare()
env.Help('LSST Association Pipeline')

