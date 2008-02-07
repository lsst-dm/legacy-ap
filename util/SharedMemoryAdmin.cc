// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Tool for inspecting and manipulating the contents of the shared
 *          memory chunk manager used by the association pipeline.
 *
 * @ingroup associate
 */

#include <cstdlib>
#include <iostream>

#include <boost/program_options.hpp>

#include <lsst/ap/Common.h>
#include <lsst/ap/Exceptions.h>
#include <lsst/ap/ChunkManager.h>


int main(int argc, char * argv[]) {

    using namespace boost::program_options;
    using namespace lsst::ap;

    try {

        options_description general("General options");
        general.add_options()
            ("help,h", "print usage help")
            ("name,n", value<std::string>()->default_value("test"),
                "the name of the shared memory object to inspect or manipulate");

        options_description inspect("Inspecting the AP chunk manager");
        inspect.add_options()
            ("visits,V", "Shows the list of visits under management")
            ("chunks,C", "Shows the list of chunks currently owned by a visit")
            ("visit,v", value<int64_t>(),
                "Displays information about the given visit, including the list of chunks owned by it")
            ("chunk,c", value<int64_t>(), "Displays information for the given chunk");

        options_description manipulate("Manipulating the AP chunk manager");
        manipulate.add_options()
            ("rollback,r", value<int64_t>(), "Rolls back the given visit [dangerous]")
            ("unlink,u", "Unlink the shared memory object underlying the AP chunk manager");

        options_description all;
        all.add(general).add(inspect).add(manipulate);

        variables_map vm;
        store(parse_command_line(argc, argv, all), vm);
        std::string name(vm["name"].as<std::string>());
        if (vm.count("help")) {
            std::cout << all;
            return EXIT_SUCCESS;
        }
        if (vm.count("visits")) {
            SharedSimpleObjectChunkManager manager(name);
            manager.printVisits(std::cout);
        }
        if (vm.count("chunks")) {
            SharedSimpleObjectChunkManager manager(name);
            manager.printChunks(std::cout);
        }
        if (vm.count("visit")) {
            SharedSimpleObjectChunkManager manager(name);
            manager.printVisit(vm["visit"].as<int64_t>(), std::cout);
        }
        if (vm.count("chunk")) {
            SharedSimpleObjectChunkManager manager(name);
            manager.printChunk(vm["chunk"].as<int64_t>(), std::cout);
        }
        if (vm.count("rollback")) {
            SharedSimpleObjectChunkManager manager(name);
            manager.endVisit(vm["rollback"].as<int64_t>(), true);
        }
        if (vm.count("unlink")) {
            SharedSimpleObjectChunkManager::destroyInstance(name);
        }

        return EXIT_SUCCESS;

    } catch (lsst::mwi::exceptions::ExceptionStack & exs) {
        std::cout << "Caught lsst::mwi::exceptions::ExceptionStack\n\t" <<
                     exs.what() << '\n' << exs.getStack()->toString("", true) << std::endl;
    } catch (std::exception & ex) {
        std::cout << "Caught std::exception\n\t" << ex.what() << std::endl;
    }
    return EXIT_FAILURE;
}

