#include <iostream>
#include <vector>

#include <libsgp4/DateTime.h>
#include <libsgp4/Eci.h>
#include <libsgp4/SGP4.h>
#include <libsgp4/Tle.h>

#include "Atom/convertCartesianStateToTwoLineElements.hpp"

typedef std::vector< double > Vec;

int main()
{
    std::cout << "Hello!" << std::endl;
    std::cout << std::endl;
    std::cout.precision( 15 );

    Vec cartesianState( 6 );
    cartesianState[ 0 ] = 3126974.99 / 1000.0;
    cartesianState[ 1 ] = -6374445.74 / 1000.0;
    cartesianState[ 2 ] = 28673.59 / 1000.0;
    cartesianState[ 3 ] = -254.91197 / 1000.0;
    cartesianState[ 4 ] = -83.30107 / 1000.0;
    cartesianState[ 5 ] = 7485.70674 / 1000.0;

    for ( unsigned int i = 0; i < 6; ++i )
    {
        std::cout << cartesianState[ i ] << std::endl;
    }
    std::cout << std::endl;

    DateTime epoch( 2014, 6, 13, 14, 59, 21 );
    std::cout << epoch << std::endl;
    std::cout << std::endl;

    Tle virtualTle = atom::convertCartesianStateToTwoLineElements< double >( cartesianState, epoch );

    std::cout << virtualTle << std::endl;
    std::cout << std::endl;

    SGP4 sgp4( virtualTle );
    Eci state = sgp4.FindPosition( 0.0 );

    std::cout << state.Position( ).x - cartesianState[ 0 ] << std::endl;
    std::cout << state.Position( ).y - cartesianState[ 1 ] << std::endl;
    std::cout << state.Position( ).z - cartesianState[ 2 ] << std::endl;
    std::cout << state.Velocity( ).x - cartesianState[ 3 ] << std::endl;
    std::cout << state.Velocity( ).y - cartesianState[ 4 ] << std::endl;
    std::cout << state.Velocity( ).z - cartesianState[ 5 ] << std::endl;
    std::cout << std::endl;

    return 0;
}