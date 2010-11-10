# -*- python -*-
import glob, os.path, re, os, sys
import lsst.SConsUtils as scons

# Custom configure tests
visCheckSrc = """
    __attribute__((visibility("hidden")))  void hiddenFunc() {}
    __attribute__((visibility("default"))) void defaultFunc() {}
    int main() {
        hiddenFunc();
        defaultFunc();
        return 0;
    }
    """

alignedCheckSrc = """
    double a[4] __attribute__((aligned(32)));
    int main() {
        return 0;
    }
    """

popcountCheckSrc = """
    int main() {
        unsigned long long ull = 0;
        unsigned long      ul  = 0;
        unsigned int       ui  = 0;
        return __builtin_popcount(ui) + __builtin_popcountl(ul) + __builtin_popcountll(ull);
    }
    """

boostInt64IsLongCheckSrc = """
    #include "boost/cstdint.hpp"
    #include "boost/static_assert.hpp"
    #include "boost/type_traits/is_same.hpp"

    int main() {
        BOOST_STATIC_ASSERT((boost::is_same<long, boost::int64_t>::value));
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

def IsGccBelow4(context):
    context.Message("Checking if CC is a version of gcc prior to 4.x ...")
    result = context.TryAction(["%s -dumpversion | grep \"^[0-3]\\..*\"" % env['CC']])[0]
    context.Result(result)
    return result

# Direct and indirect dependencies of ap
dependencies = ["boost",
                "python",
                "mysqlclient",
                "wcslib",
                "minuit2",
                "base",
                "pex_exceptions",
                "utils",
                "daf_base",
                "pex_logging",
                "security",
                "pex_policy",
                "daf_persistence",
                "daf_data",
                "afw",
                "mops"]

#
# Setup our environment
#
env = scons.makeEnv("ap",
                    r"$HeadURL$",
                    [["boost", "boost/version.hpp", "boost_system:C++"],
                     ["boost", "boost/filesystem.hpp", "boost_filesystem:C++"],
                     ["boost", "boost/regex.hpp", "boost_regex:C++"],
                     ["boost", "boost/serialization/base_object.hpp", "boost_serialization:C++"],
                     ["boost", "boost/program_options.hpp", "boost_program_options:C++"],                    
                     ["boost", "boost/test/unit_test.hpp", "boost_unit_test_framework:C++"],
                     ["python", "Python.h"],
                     ["mysqlclient", "mysql/mysql.h", "mysqlclient:C++"],
                     ["wcslib", "wcslib/wcs.h", "m wcs"],
                     ["minuit2", "Minuit2/FCNBase.h", "Minuit2:C++"],
                     ["gsl", "gsl/gsl_rng.h", "gslcblas gsl"],
                     ["eigen", "Eigen/Core.h"],
                     ["base", "lsst/base.h"],
                     ["pex_exceptions", "lsst/pex/exceptions.h", "pex_exceptions:C++"],
                     ["utils", "lsst/utils/Utils.h", "utils:C++"],
                     ["daf_base", "lsst/daf/base.h", "daf_base:C++"],
                     ["pex_logging", "lsst/pex/logging/Trace.h", "pex_logging:C++"],
                     ["security", "lsst/security/Security.h", "security:C++"],
                     ["pex_policy", "lsst/pex/policy/Policy.h", "pex_policy:C++"],
                     ["daf_persistence", "lsst/daf/persistence.h", "daf_persistence:C++"],
                     ["daf_data", "lsst/daf/data.h", "daf_data:C++"],
                     ["afw", "lsst/afw/detection/DiaSource.h", "afw:C++"],
                     ["xpa", "xpa.h", "xpa", "XPAPuts"],
                     ["mops", "lsst/mops/MovingObjectPrediction.h"]
                    ])

env.Help("""
LSST Association Pipeline
""")

# Libraries needed to link libraries/executables
pkg = env["eups_product"]
env.libs[pkg] += env.getlibs(" ".join(dependencies))

sysLibs = ['z', 'pthread']
if env['PLATFORM'] == 'posix':
    sysLibs += ['rt']
env.libs[pkg] += sysLibs

# Run configure tests
if not env.CleanFlagIsSet():
    conf = Configure(env, custom_tests = {'CustomCompilerFlag' : CustomCompilerFlag,
                                          'CustomCompileCheck' : CustomCompileCheck,
                                          'CustomLinkCheck'    : CustomLinkCheck,
                                          'IsGccBelow4'        : IsGccBelow4})
    if env['PLATFORM'] == 'posix':
        # POSIX platforms have AIO functionality in librt
        if not conf.CheckLibWithHeader('rt', 'aio.h', 'C'):
            print 'Missing support for Posix AIO'
            Exit(1)

    # If one of the randomized unit tests fail, uncomment and
    # set the following defines to obtain repeatable behaviour
    #conf.env.Append(CPPFLAGS = ' -DLSST_AP_INIT_RANDOM_SEC=')
    #conf.env.Append(CPPFLAGS = ' -DLSST_AP_INIT_RANDOM_NSEC=')

    # Required for [U]INT64_C()
    conf.env.Append(CPPFLAGS = ' -D__STDC_CONSTANT_MACROS')
    # Indicate that shared libraries are being built
    conf.env.Append(CPPFLAGS = ' -DLSST_AP_SHARED_LIBRARY_BUILD=1')
    # compiler flags and features
    if not conf.IsGccBelow4():
        conf.CustomCompilerFlag('-fvisibility-inlines-hidden')
    if conf.CustomCompileCheck('Checking for __attribute__((visibility)) support... ', visCheckSrc):
        conf.env.Append(CPPFLAGS = ' -DLSST_AP_HAVE_VISIBILITY=1')
    if conf.CustomCompileCheck('Checking for __attribute__((aligned)) support... ', alignedCheckSrc):
        conf.env.Append(CPPFLAGS = ' -DLSST_AP_HAVE_ALIGNED=1')
    if conf.CustomCompileCheck('Checking for __builtin_popcount... ', popcountCheckSrc):
        conf.env.Append(CPPFLAGS = ' -DLSST_AP_HAVE_BUILTIN_POPCOUNT=1')
    # Without some help, SWIG disagrees with boost on the actual type of int64_t
    if conf.CustomCompileCheck('Checking whether boost::int64_t is long ... ',
                               boostInt64IsLongCheckSrc, extension='.cc'):
        conf.env.Append(SWIGFLAGS = '-DSWIGWORDSIZE64')
    # Platform features
    if conf.CheckFunc('clock_gettime'): # Linux/Solaris: prototype in <time.h>
        conf.env.Append(CPPFLAGS = ' -DLSST_AP_HAVE_CLOCK_GETTIME=1')
    env = conf.Finish()


# Build/install things
directories = ["bin",
               "lib",
               "python/lsst/" + pkg + "/cluster",
               "python/lsst/" + pkg + "/match",
               "tests",
               "doc"]

for d in directories:
    try:
        SConscript(os.path.join(d, "SConscript"))
    except Exception, e:
        print >> sys.stderr, "%s: %s" % (os.path.join(d, "SConscript"), e)

env['IgnoreFiles'] = r"(~$|\.pyc$|^\.svn$|\.o$)"

Alias("install", [env.Install(env['prefix'], "python"),
                  env.Install(env['prefix'], "include"),
                  env.Install(env['prefix'], "lib"),
                  env.Install(env['prefix'], "policy"),
                  env.Install(env['prefix'], "bin"),
                  env.Install(env['prefix'], "doc"),
                  env.InstallEups(env['prefix'] + "/ups")])

scons.CleanTree(r"*~ core *.so *.os *.o")

# Build TAGS files
files = scons.filesToTag()
if files:
    env.Command("TAGS", files, "etags -o $TARGET $SOURCES")

env.Declare()

