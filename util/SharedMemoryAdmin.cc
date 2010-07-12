// -*- lsst-c++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 

/**
 * @file
 * @brief   Tool for inspecting and manipulating the contents of the shared
 *          memory chunk manager used by the association pipeline.
 *
 * @ingroup associate
 */

#include <cstdlib>
#include <iostream>

#include "boost/program_options.hpp"

#include "lsst/ap/Common.h"
#include "lsst/ap/ChunkManager.h"


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
            ("visit,v", value<int>(),
                "Displays information about the given visit, including the list of chunks owned by it")
            ("chunk,c", value<int>(), "Displays information for the given chunk");

        options_description manipulate("Manipulating the AP chunk manager");
        manipulate.add_options()
            ("rollback,r", value<int>(), "Rolls back the given visit [dangerous]")
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
            SharedObjectChunkManager manager(name);
            manager.printVisits(std::cout);
        }
        if (vm.count("chunks")) {
            SharedObjectChunkManager manager(name);
            manager.printChunks(std::cout);
        }
        if (vm.count("visit")) {
            SharedObjectChunkManager manager(name);
            manager.printVisit(vm["visit"].as<int>(), std::cout);
        }
        if (vm.count("chunk")) {
            SharedObjectChunkManager manager(name);
            manager.printChunk(vm["chunk"].as<int>(), std::cout);
        }
        if (vm.count("rollback")) {
            SharedObjectChunkManager manager(name);
            manager.endVisit(vm["rollback"].as<int>(), true);
        }
        if (vm.count("unlink")) {
            SharedObjectChunkManager::destroyInstance(name);
        }

        return 0;

    } catch (lsst::pex::exceptions::Exception & ex) {
        std::cout << "Caught lsst::pex::exceptions::Exception :\n\n" << ex;
    } catch (std::exception & ex) {
        std::cout << "Caught std::exception : " << ex.what() << std::endl;
    }
    return 1;
}

