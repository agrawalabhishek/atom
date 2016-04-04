// Program to diagnose the cartesian state to two line element converter.
// Copyright (c) 2016, Abhishek Agrawal (abhishek.agrawal@protonmail.com)

#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include <libsgp4/DateTime.h>
#include <libsgp4/Eci.h>
#include <libsgp4/Globals.h>
#include <libsgp4/SGP4.h>
#include <libsgp4/Tle.h>
#include <libsgp4/Vector.h>

#include <Astro/astro.hpp>

#include <boost/array.hpp>

#include <SML/sml.hpp>

#include "Atom/convertCartesianStateToTwoLineElements.hpp"

typedef boost::array< double, 3 > Vector3;
typedef boost::array< double, 6 > Vector6;

int main( void )
{
    std::cout << "Execute Cartesian State To Two-Line Element Converter Diagnostics..."
    std::cout.precision( 15 );

    // Generate random keplerian elements


    return 0;
}
