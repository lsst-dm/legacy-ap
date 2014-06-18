// -*- lsst-c++ -*-

/*
 * LSST Data Management System
 * Copyright 2012 LSST Corporation.
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

/** @file
  * @brief  SourceProcessingControl implementation.
  * @author Serge Monkewitz
  */
#include "lsst/ap/cluster/SourceProcessingControl.h"

#include "lsst/pex/exceptions.h"

using lsst::pex::exceptions::InvalidParameterError;


namespace lsst { namespace ap { namespace cluster {

SourceProcessingControl::SourceProcessingControl() :
    exposurePrefix("exposure"),
    clusterPrefix("cluster"),
    multiBand(false),
    coadd(false),
    badFlagFields(),
    // The default flux scale S is chosen such that the AB magnitude
    // -2.5*log10(S*F/F_0) - 48.6 = 0 when F = F_0, where F is an
    // uncalibrated flux and F_0 is the flux of a 0-magnitude source.
    fluxScale(3.63078054770101342467371212362e-20),
    fluxUnit("erg/s/cm^2/Hz"),
    fluxFields(),
    shapeFields()
{
    fluxFields.push_back("flux.gaussian");
    fluxFields.push_back("flux.naive");
    fluxFields.push_back("flux.psf");
    fluxFields.push_back("flux.sinc");
    fluxFields.push_back("flux.kron");
    fluxFields.push_back("multishapelet.exp.flux");
    fluxFields.push_back("multishapelet.dev.flux");
    fluxFields.push_back("multishapelet.combo.flux");
    shapeFields.push_back("shape.sdss");
}

SourceProcessingControl::~SourceProcessingControl() { }

}}} // namespace lsst::ap::cluster
